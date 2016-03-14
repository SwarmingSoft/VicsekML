//entropy.cpp v17 05.02.2016

//EntropyVicsek.exe data entropy staticentropy 256 16. 0.01 topological 1 20 1 9999 109999 100

//CARE!!!: eigen's .inverse() does not work with "-o3"
//confer: http://eigen.tuxfamily.org/bz/show_bug.cgi?id=424

#include <string>
#include <sstream>
#include <array>
#include <iostream>
#include <cmath>

#include <dlib/optimization.h>
#include <dlib/matrix.h>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

#include "constants.hpp"
#include "tests.hpp"
#include "MyMatrix.hpp"

enum class Model {DynamicEntropy, StaticEntropy, DynamicOptimization, GeneralizedEntropy, GeneralizedEntropyJ, GeneralizedEntropySingle};

class CMDParameterEntropy
{

public:

    std::string in_file, out_file;
    Model eval_model;
    unsigned int sim_N;
    real sim_L;
    real sim_deltat;

    Geometry eval_geometry;
    real param_start, param_stop, param_step;

    unsigned int time_start, time_stop, time_step;

    CMDParameterEntropy(int argc, char** argv)
    {
        if (argc == 1)
        {
            std::cerr << "inputfile(string)\toutputfile(string)\tmodel(dynamicentropy,staticentropy,dynamicoptimization,generalizedentropy,generalizedentropyJ)\tN(int)\tL(real)\tdeltat(real)\tGeometry(metric,topological,voronoi)\tparam_start(real)\tparam_stop(real)\tparam_step(real)\ttime_start(int)\ttime_stop(int)\ttime_step(int)" << std::endl;
            throw std::runtime_error("Parameters missing");
        }
        assert(argc == 1 + 13);

        in_file = std::string(argv[1]) + ".txt";
        out_file = std::string(argv[2]) + ".txt";

        std::string _model = argv[3];
        if (_model == "dynamicentropy")
        {
            eval_model = Model::DynamicEntropy;
        }
        else if (_model == "staticentropy")
        {
            eval_model = Model::StaticEntropy;
        }
        else if (_model == "dynamicoptimization")
        {
            eval_model = Model::DynamicOptimization;
        }
        else if (_model == "generalizedentropy")
        {
            eval_model = Model::GeneralizedEntropy;
        }
        else if (_model == "generalizedentropyJ")
        {
            eval_model = Model::GeneralizedEntropyJ;
        }
        else
        {
            throw std::runtime_error("Invalid model parameter");
        }

        sim_N = std::atoi(argv[4]);
        sim_L = std::atof(argv[5]);
        sim_deltat = std::atof(argv[6]);

        std::string _geometry = argv[7];

        if (_geometry == "topological")
        {
            eval_geometry = Geometry::topological;
            param_start = std::atof(argv[8]);
            param_stop = std::atof(argv[9]);
            param_step = std::atof(argv[10]);
        }
        else if (_geometry == "metric")
        {
            eval_geometry = Geometry::metric;
            param_start = std::atof(argv[8]);
            param_stop = std::atof(argv[9]);
            param_step = std::atof(argv[10]);
        }
        else if (_geometry == "voronoi")
        {
            eval_geometry = Geometry::voronoi;
            param_start = 6;
            param_stop = 7;
            param_step = 1;
        }
        else
        {
            throw std::runtime_error("Invalid geometry parameter");
        }

        time_start = std::atoi(argv[11]);
        time_stop = std::atoi(argv[12]);
        time_step = std::atoi(argv[13]);
    }

    std::string header(void)
    {
        std::stringstream ss;

        ss << "In-file: " << in_file << " Model: " << int(eval_model) << " N: " << sim_N << " L: " << sim_L << " deltat: " << sim_deltat << " Geometry: " << int(eval_geometry) << " Param-start: " << param_start << " Param-stop: " << param_stop << " Param-step: " << param_step << " Time-start: " << time_start << " Time-stop: " << time_stop << " Time-step: " << time_step;
        return ss.str();
    }

};

// Dlib helper

typedef dlib::matrix<real,0,1> dlib_vector;

std::string StringVector(const dlib_vector &v)
{
    std::stringstream ss;
    for (int i = 0; i < v.size() - 1; ++i)
    {
        ss << v(i) << " ";
    }
    ss << v(v.size()-1);
    return ss.str();
}

// Eigen helper

typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixReal;
typedef Eigen::Matrix<real, Eigen::Dynamic, 1> EigenVectorReal;

void ConvertEigenMatrixToRealMatrix(const EigenMatrixReal &a, RealMatrix &b)
{
    assert(a.rows() == b.GetRows() && a.cols() == b.GetColumns());
    unsigned int r, c;

    for (r = 0; r < a.rows(); ++r)
    {
        for (c = 0; c < a.cols(); ++c)
        {
            b(r, c) = a(r, c);
        }
    }
}

// RealMatrix helper

RealMatrix convert_to_real_colvector(const std::vector<real> &v)
{
    RealMatrix c(v.size(),1);

    for (unsigned int i = 0; i < v.size(); ++i)
    {
        c(i, 0) = v[i];
    }

    return c;
}

VectorMatrix convert_to_vector_colvector(const std::vector<Vector> &v)
{
    VectorMatrix c(v.size(),1);

    for (unsigned int i = 0; i < v.size(); ++i)
    {
        c(i, 0) = v[i];
    }

    return c;
}

//Optimization

class BufferedN
{
    real L;
    std::map<real, RealMatrix> ns;
    std::map<real, RealMatrix> nis;
    real paramstart, paramend, paramstep;
    const std::vector<Vector> &positions;
    Geometry geometry;

public:

    const unsigned int N;

    BufferedN(real L, const std::vector<Vector> &positions) : L(L), positions(positions), N(positions.size())
    {
        //pass
    }

    void InitTopological(real paramstart, real paramend, real paramstep)
    {
        geometry = Geometry::topological;
        this->paramstart = paramstart;
        this->paramend = paramend;
        this->paramstep = paramstep;
    }

    void InitMetric(real paramstart, real paramend, real paramstep)
    {
        geometry = Geometry::metric;
        this->paramstart = paramstart;
        this->paramend = paramend;
        this->paramstep = paramstep;
    }

    void InitVoronoi(void)
    {
        geometry = Geometry::voronoi;
        paramstart = 6;
        paramend = 7;
        paramstep = 1;
    }

    void Calculate(void)
    {
        MATRIX n(N, N);
        RealMatrix m(N, N);

        std::vector<real> ni(N);
        RealMatrix mi(N,1);
        for (real param = paramstart; param < paramend; param += paramstep)
        {
            switch (geometry)
            {
                case Geometry::voronoi:
                    calculate_n_voronoi(n, positions, L);
                    assert(is_symmetric(n));
                    break;
                case Geometry::metric:
                    calculate_n_metric(n, positions, L, param);
                    assert(is_symmetric(n));
                    break;
                case Geometry::topological:
                    calculate_n_topological(n, positions, L, param);
                    assert(is_number_of_neighbours_equal(n));
                    break;
            }

            create_ni(n, ni);

            // std::vector<real> ni => RealMatrix mi
            nis[param] = convert_to_real_colvector(ni);

            // MATRIX n => RealMatrix m
            m = n;
            ns[param] = m;
        }
    }

    real empnc(real param) const
    {
        const RealMatrix &v = nis.at(param);

        assert(v.GetColumns() == 1);
        assert(v.GetRows() > 0);

        real sum = 0;
        for (unsigned int i = 0; i < v.GetRows(); ++i)
        {
            sum += v(i,0);
        }
        return sum / v.GetRows();
    }

    RealMatrix ij(real param) const
    {
        /*unsigned int discrete_neighbours = std::round(param);

        if (discrete_neighbours >= start && discrete_neighbours < end)
        {
            return ns.at(discrete_neighbours - start);
        }
        else
        {
            //return dlib::zeros_matrix<real>(N,N);
            return RealMatrix::Zero(N,N);
        }*/
        //std::cout << "Read ns: " << param << std::endl;
        return ns.at(param);
    }

    RealMatrix i(real param) const
    {
        /*unsigned int discrete_neighbours = std::round(param);

        if (discrete_neighbours >= start && discrete_neighbours < end)
        {
            return nis.at(discrete_neighbours - start);
        }
        else
        {
            //return dlib::zeros_matrix<real>(N,1);
            return RealMatrix::Zero(N,1);
        }*/
        //std::cout << "Read nis: " << param << std::endl;
        return nis.at(param);
    }
};

class DynamicLikelihood
{
    real deltat;
    BufferedN n;
    VectorMatrix pi, pip;
    real param = 0;

    unsigned int N;

public:

    DynamicLikelihood(real deltat, BufferedN n, VectorMatrix pi, VectorMatrix pip) : deltat(deltat), n(n), pi(pi), pip(pip), N(n.N)
    {
        //pass
    }

    //pi'*pi' = sum_i pi'_i^2
    /*real sum1(const VectorMatrix &pip) const
    {
        int i;
        real sum = 0;
        real tmp;

        for (i = 0; i < pip.size(); ++i)
        {
            tmp = pip(i)*pip(i);
            assert(-1. <= tmp && tmp <= 1.);
            sum += tmp;
        }

        return sum;
    }

    //pi^T M^T pi' == pi' M pi = sum_ij M_ij pi'_i pi_j
    real sum2(const RealMatrix &M, const VectorMatrix &pi, const VectorMatrix &pip) const
    {
        int i, j;
        real sum = 0;
        real tmp;

        for (i = 0; i < pi.size(); ++i)
        {
            for (j = 0; j < pi.size(); ++j)
            {
                assert(0 <= i && i < N && 0 <= j && j < N);
                tmp = pip(i) * pi(j);
                //assert(-real(N) <= M(i,j) && M(i,j) <= real(N));
                if (!(-real(N) <= M(i,j) && M(i,j) <= real(N)))
                {
                    std::cout << "!" << -real(N) << "<=" << M(i,j) << "<=" << real(N);
                }
                assert(-1. <= tmp && tmp <= 1.);
                sum += M(i,j) * tmp;
            }
        }

        return sum;
    }


    // pi^T M^T M pi = sum_ijk M_ij pi_j M_ik pi_k
    real sum3(const RealMatrix &M, const VectorMatrix &pi) const
    {
        int i, j, k;
        real sum = 0;
        real tmp;

        for (i = 0; i < pi.size(); ++i)
        {
            for (j = 0; j < pi.size(); ++j)
            {
                for (k = 0; k < pi.size(); ++k)
                {
                    tmp = pi(j) * pi(k);
                    assert(-real(N) <= M(i,j) && M(i,j) <= real(N));
                    assert(-real(N) <= M(i,k) && M(i,k) <= real(N));
                    assert(-1. <= tmp && tmp <= 1.);
                    sum += M(i,j) * M(i,k) * tmp;
                }
            }
        }

        return sum;
    }*/

    /* sum version:
        real s1 = sum1(pip);
        real s2 = sum2(M, pi, pip);
        real s3 = sum3(M, pi);

        real p1 = - std::log(2*deltat*T*2*PI);
        real p2 = (s1 -2*s2 + s3)/N;

        return -(d-1.)/2. * p1 + 0.5 * 1./(2*deltat*T)*p2;
    */

    void set_param(real param)
    {
        this->param = param;
    }

    real operator()(const dlib_vector &m) const
    {
        // -N STAYS UNSIGNED
        const real J = m(0);
        const real sigma2 = m(1);
        //const real n_c = m(2);

        RealMatrix A = 1./(deltat*sigma2)*RealMatrix::Identity(N);
        RealMatrix M = RealMatrix::Identity(N) - J*deltat*RealMatrix::Diag(n.i(param)) + J*deltat*n.ij(param);

        real ret = -(N/2.) * std::log(2*PI*deltat*sigma2) - real((pip-M*pi).Transposed()*A*(pip-M*pi)) / 2.;

        //std::cout << "1: " << -(N/2.) * std::log(2*PI*deltat*sigma2) << " 2: " << - real((pip-M*pi).Transposed()*A*(pip-M*pi)) / 2. << std::endl;
        //std::cout << "J: " << J << " sigma2: " << sigma2 << " param: " << param << " => " << ret << std::endl;

        return -ret; //invert for minimization
    }

    real operator__(const dlib_vector &m) const
    {
        return m(0)+m(1)+m(2);
    }
};

// general optimization

class BufferedLambda
{
    real L;
    std::map<real, EigenMatrixReal> ns;
    std::map<real, RealMatrix> nis;
    real paramstart, paramend, paramstep;
    const std::vector<Vector> &positions;
    Geometry geometry;

public:

    const unsigned int N;

    BufferedLambda(real L, const std::vector<Vector> &positions) : L(L), positions(positions), N(positions.size())
    {
        //pass
    }

    void InitTopological(real paramstart, real paramend, real paramstep)
    {
        geometry = Geometry::topological;
        this->paramstart = paramstart;
        this->paramend = paramend;
        this->paramstep = paramstep;
    }

    void InitMetric(real paramstart, real paramend, real paramstep)
    {
        geometry = Geometry::metric;
        this->paramstart = paramstart;
        this->paramend = paramend;
        this->paramstep = paramstep;
    }

    void InitVoronoi(void)
    {
        geometry = Geometry::voronoi;
        paramstart = 6;
        paramend = 7;
        paramstep = 1;
    }

    void Calculate(void)
    {
        MATRIX n(N, N);
        RealMatrix m(N, N);

        std::vector<real> ni(N);
        RealMatrix mi(N,1);
        for (real param = paramstart; param < paramend; param += paramstep)
        {
            switch (geometry)
            {
                case Geometry::voronoi:
                    calculate_n_voronoi(n, positions, L);
                    assert(is_symmetric(n));
                    break;
                case Geometry::metric:
                    calculate_n_metric(n, positions, L, param);
                    assert(is_symmetric(n));
                    break;
                case Geometry::topological:
                    calculate_n_topological(n, positions, L, param);
                    assert(is_number_of_neighbours_equal(n));
                    break;
            }

            create_ni(n, ni);

            EigenMatrixReal lambda(N,N);

            unsigned int r, c;
            for (r = 0; r < n.GetHeight(); ++r)
            {
                for (c = 0; c < n.GetWidth(); ++c)
                {
                    if (r == c)
                    {
                        //n_i delta_ij - n_ij
                        lambda(r, c) = ni[r] - *n.Get(r, c);
                    }
                    else
                    {
                        lambda(r, c) = - *n.Get(r, c);
                    }
                }
            }

            assert(is_weakly_diagonally_dominant_and_non_negativ_diagonal_entries(lambda));

            // std::vector<real> ni => RealMatrix mi
            nis[param] = convert_to_real_colvector(ni);

            // MATRIX n => RealMatrix m
            ns[param] = lambda;
        }
    }

    real empnc(real param) const
    {
        const RealMatrix &v = nis.at(param);

        assert(v.GetColumns() == 1);
        assert(v.GetRows() > 0);

        real sum = 0;
        for (unsigned int i = 0; i < v.GetRows(); ++i)
        {
            sum += v(i,0);
        }
        return sum / v.GetRows();
    }

    EigenMatrixReal ij(real param) const
    {
        return ns.at(param);
    }

    RealMatrix i(real param) const
    {
        return nis.at(param);
    }
};

real LogTrace(const EigenMatrixReal& m)
{
    assert(m.rows() == m.cols());
    real ret = 0;
    for (unsigned int i = 0; i < m.rows(); ++i)
    {
        ret += std::log(m(i, i));
    }
    return ret;
}


class GeneralLikelihoodSingle
{
    real deltat;
    BufferedLambda Blambda;
    VectorMatrix pi, pip;
    real param = 0;
    unsigned int N;

public:

    GeneralLikelihoodSingle(real deltat, BufferedLambda Blambda, VectorMatrix pi, VectorMatrix pip) : deltat(deltat), Blambda(Blambda), pi(pi), pip(pip), N(Blambda.N)
    {

    }

    void set_param(real param)
    {
        this->param = param;
    }

    real operator()(const dlib_vector &m) const
    {
        const real J = m(0);
        const real sigma2 = m(1);
        const EigenMatrixReal &lambda = Blambda.ij(param);

        assert(lambda == lambda.transpose());
        Eigen::SelfAdjointEigenSolver<EigenMatrixReal> es;
        es.compute(lambda, /* computeEigenvectors = */ Eigen::ComputeEigenvectors);
        assert(es.info() == Eigen::ComputationInfo::Success);

        auto eval_vector = es.eigenvalues();
        //std::cout << eval_vector << std::endl;
        EigenMatrixReal trafo = es.eigenvectors(); // consists of eigen vectors


        //create intexpdiag
        EigenMatrixReal intexpdiag = EigenMatrixReal::Zero(N, N); //set all entries to zero
        for (unsigned int i = 0; i < N; ++i)
        {
            assert(eval_vector(i) > -0.0001);
            if (eval_vector(i) < 0.0001)
            {
                intexpdiag(i, i) = sigma2 * deltat;
            }
            else
            {
                intexpdiag(i, i) = sigma2 *  (1. - std::exp(-2*J*eval_vector(i)*deltat)) / (2 * eval_vector(i) * J);
            }
        }

        //create inverse of intexpdiag
        EigenMatrixReal intexpdiag_inv = EigenMatrixReal::Zero(N, N); //set all entries to zero
        for (unsigned int i = 0; i < N; ++i)
        {
            assert(intexpdiag(i, i) != real(0));
            intexpdiag_inv(i, i) = 1. / intexpdiag(i, i);
        }

        EigenMatrixReal xiinv_eig = trafo * (intexpdiag_inv * trafo.inverse());

        RealMatrix xiinv(N, N);
        ConvertEigenMatrixToRealMatrix(xiinv_eig, xiinv);

        // Tr(X) == Tr(U X U^-1) see wiki
        real logdet = LogTrace(intexpdiag);

        EigenMatrixReal propagator_eig = (-J * deltat * lambda).exp();

        //RealMatrix lambda_our(N, N);
        //ConvertEigenMatrixToRealMatrix(lambda, lambda_our);
        //RealMatrix propagator = MatrixExponential(-J * deltat * lambda_our, 4);
        //std::cout << propagator_eig << std::endl;
        //PrintMatrix(propagator);

        RealMatrix propagator(N, N);
        ConvertEigenMatrixToRealMatrix(propagator_eig, propagator);

        real ll_logt = -(N/2.) * std::log(2*PI);
        real ll_logdet = -logdet / 2.;
        real ll_mpi = -real((pip-propagator*pi).Transposed()*(xiinv*(pip-propagator*pi))) / 2.;

        real ll = ll_logt + ll_logdet + ll_mpi;
        //std::cout << "Likelihood: " << ll << " J: " << J << " T: " << T << " logt: " << ll_logt << " logdet: " << ll_logdet << " mpi: " << ll_mpi << std::endl;
        return -ll; //invert for minimization
    }
};

class GeneralLikelihoodSum
{
protected:
    real deltat;
    std::vector<unsigned int> times;
    Pi pis;
    unsigned int N;

    std::vector<EigenMatrixReal> Lambda;

    std::vector<EigenMatrixReal> U;
    std::vector<EigenMatrixReal> UI;
    std::vector<EigenMatrixReal> ev; // EigenMatrixReal?? //typedef Matrix< double , Dynamic , 1> VectorXd

    std::vector<std::vector<real>> nis; // only for w/o J approx

public:

    std::vector<real> empncs;

    GeneralLikelihoodSum(real deltat, const std::vector<unsigned int> &times, const Pi &pis) : deltat(deltat), times(times), pis(pis), N(pis.Get(times[0]).size()) //pis.size() should work too..
    {
        empncs.resize(times.size());
        Lambda.resize(times.size());
        ev.resize(times.size());
        U.resize(times.size());
        UI.resize(times.size());
        nis.resize(times.size()); // only for w/o J approx
        //std::cout << pis.Get(times[0]).size() << std::endl;
    }

    void Init(const Posdict &positions, real param, Geometry geometry, real L)
    {
        //reserve for all vectors?????
        //pis.size() == len(times) and not == N ????? --> pis.Get(times[0]).size()

        //empncs.resize(0);
        //Lambda.resize(0);
        //ev.resize(0);
        //U.resize(0);
        //UI.resize(0);
        for (unsigned int t = 0; t < times.size(); ++t)
        {
            MATRIX n(N, N);
            RealMatrix m(N, N);

            std::vector<real> ni(N);
            RealMatrix mi(N,1);

            switch (geometry)
            {
                case Geometry::voronoi:
                    calculate_n_voronoi(n, positions.at(times[t]), L);
                    assert(is_symmetric(n));
                    break;
                case Geometry::metric:
                    calculate_n_metric(n, positions.at(times[t]), L, param);
                    assert(is_symmetric(n));
                    break;
                case Geometry::topological:
                    calculate_n_topological(n, positions.at(times[t]), L, param);
                    assert(is_number_of_neighbours_equal(n));
                    break;
            }

            create_ni(n, ni);

            EigenMatrixReal lambda(N,N);

            unsigned int r, c;
            for (r = 0; r < n.GetHeight(); ++r)
            {
                for (c = 0; c < n.GetWidth(); ++c)
                {
                    if (r == c)
                    {
                        //n_i delta_ij - n_ij
                        lambda(r, c) = ni[r] - *n.Get(r, c);
                    }
                    else
                    {
                        lambda(r, c) = - *n.Get(r, c);
                    }
                }
            }

            assert(is_weakly_diagonally_dominant_and_non_negativ_diagonal_entries(lambda));

            empncs[t] = average(ni);
            //empncs.push_back(average(ni));
            Lambda[t] = lambda;
            //Lambda.push_back(lambda);
            nis[t] = ni; // only for w/o J approx
        }
    }

    void CalculateUUI(void)
    {
        Eigen::SelfAdjointEigenSolver<EigenMatrixReal> es;
        for (unsigned int t = 0; t < times.size(); ++t)
        {
            es.compute(Lambda[t], /* computeEigenvectors = */ Eigen::ComputeEigenvectors);
            assert(es.info() == Eigen::ComputationInfo::Success);

            ev[t] = es.eigenvalues();
            //ev.push_back(es.eigenvalues());
            U[t] = es.eigenvectors(); // consists of eigen vectors
            //U.push_back(es.eigenvectors()); // consists of eigen vectors
            UI[t] = U[t].inverse();
            //UI.push_back(U[t].inverse());
        }
    }
};

class GeneralLikelihoodSumWithJApprox : public GeneralLikelihoodSum
{
public:

    GeneralLikelihoodSumWithJApprox(real deltat, const std::vector<unsigned int> &times, const Pi &pis) : GeneralLikelihoodSum(deltat, times, pis) { }

    real operator()(const dlib_vector &m) const
    {
        const real J = m(0);
        const real sigma2 = m(1);

        real ll_sum = 0.;
        for (unsigned int t = 0; t < times.size(); ++t)
        {
            VectorMatrix pi = convert_to_vector_colvector(pis.Get(times[t]));
            VectorMatrix pip = convert_to_vector_colvector(pis.Get(times[t]+1));

            //create intexpdiag
            EigenMatrixReal intexpdiag = EigenMatrixReal::Zero(N, N); //set all entries to zero
            for (unsigned int i = 0; i < N; ++i)
            {
                assert(ev[t](i) > -0.0001);
                if (ev[t](i) < 0.0001)
                {
                    intexpdiag(i, i) = sigma2 * deltat;
                }
                else
                {
                    intexpdiag(i, i) = sigma2 *  (1. - std::exp(-2.*J*ev[t](i)*deltat)) / (2. * ev[t](i) * J);
                }
            }

            //create inverse of intexpdiag
            EigenMatrixReal intexpdiag_inv = EigenMatrixReal::Zero(N, N); //set all entries to zero
            for (unsigned int i = 0; i < N; ++i)
            {
                assert(intexpdiag(i, i) != real(0));
                intexpdiag_inv(i, i) = 1. / intexpdiag(i, i);
            }

            EigenMatrixReal xiinv_eig = U[t] * (intexpdiag_inv * UI[t]);

            RealMatrix xiinv(N, N);
            ConvertEigenMatrixToRealMatrix(xiinv_eig, xiinv);

            // Tr(X) == Tr(U X U^-1) see wiki
            real logdet = LogTrace(intexpdiag);

            EigenMatrixReal propagator_eig = (-J * deltat * Lambda[t]).exp();

            //RealMatrix lambda_our(N, N);
            //ConvertEigenMatrixToRealMatrix(lambda, lambda_our);
            //RealMatrix propagator = MatrixExponential(-J * deltat * lambda_our, 4);
            //std::cout << propagator_eig << std::endl;
            //PrintMatrix(propagator);

            RealMatrix propagator(N, N);
            ConvertEigenMatrixToRealMatrix(propagator_eig, propagator);

            real ll_logt = -(N/2.) * std::log(2*PI);
            real ll_logdet = -logdet / 2.;
            real ll_mpi = -real((pip-propagator*pi).Transposed()*(xiinv*(pip-propagator*pi))) / 2.;

            real ll = ll_logt + ll_logdet + ll_mpi;
            ll_sum += ll;
        }
        //std::cout << "Likelihood: " << ll << " J: " << J << " T: " << T << " logt: " << ll_logt << " logdet: " << ll_logdet << " mpi: " << ll_mpi << std::endl;
        return -ll_sum; //invert for minimization
    }
};

class GeneralLikelihoodSumWithoutJApprox : public GeneralLikelihoodSum
{
public:

    GeneralLikelihoodSumWithoutJApprox(real deltat, const std::vector<unsigned int> &times, const Pi &pis) : GeneralLikelihoodSum(deltat, times, pis) { }

    real operator()(const dlib_vector &m) const
    {
        const real JV = m(0);
        const real sigma2 = m(1);

        real ll_sum = 0.;
        for (unsigned int t = 0; t < times.size(); ++t)
        {
            VectorMatrix pi = convert_to_vector_colvector(pis.Get(times[t]));
            VectorMatrix pip = convert_to_vector_colvector(pis.Get(times[t]+1));

            assert(N == nis[t].size());

            //create intexpdiag and J_i
            EigenMatrixReal intexpdiag = EigenMatrixReal::Zero(N, N); //set all entries to zero
            EigenVectorReal J_i(N); // uninitialized
            for (unsigned int i = 0; i < N; ++i)
            {
                J_i[i] = JV / (1 + JV*deltat*nis[t][i]);
                assert(ev[t](i) > -0.0001);
                if (ev[t](i) < 0.0001)
                {
                    intexpdiag(i, i) = sigma2 * deltat;
                }
                else
                {
                    intexpdiag(i, i) = sigma2 *  (1. - std::exp(-2.*J_i[i]*ev[t](i)*deltat)) / (2. * ev[t](i) * J_i[i]); // just use J_i instead of J
                }
            }

            //create inverse of intexpdiag
            EigenMatrixReal intexpdiag_inv = EigenMatrixReal::Zero(N, N); //set all entries to zero
            for (unsigned int i = 0; i < N; ++i)
            {
                assert(intexpdiag(i, i) != real(0));
                intexpdiag_inv(i, i) = 1. / intexpdiag(i, i);
            }

            EigenMatrixReal xiinv_eig = U[t] * (intexpdiag_inv * UI[t]);

            RealMatrix xiinv(N, N);
            ConvertEigenMatrixToRealMatrix(xiinv_eig, xiinv);

            // Tr(X) == Tr(U X U^-1) see wiki
            real logdet = LogTrace(intexpdiag);

            EigenMatrixReal propagator_eig = (-deltat * J_i.asDiagonal() * Lambda[t]).exp(); // use J*Id instead of J here

            //RealMatrix lambda_our(N, N);
            //ConvertEigenMatrixToRealMatrix(lambda, lambda_our);
            //RealMatrix propagator = MatrixExponential(-J * deltat * lambda_our, 4);
            //std::cout << propagator_eig << std::endl;
            //PrintMatrix(propagator);

            RealMatrix propagator(N, N);
            ConvertEigenMatrixToRealMatrix(propagator_eig, propagator);

            real ll_logt = -(N/2.) * std::log(2*PI);
            real ll_logdet = -logdet / 2.;
            real ll_mpi = -real((pip-propagator*pi).Transposed()*(xiinv*(pip-propagator*pi))) / 2.;

            real ll = ll_logt + ll_logdet + ll_mpi;
            ll_sum += ll;
        }
        //std::cout << "Likelihood: " << ll << " J: " << J << " T: " << T << " logt: " << ll_logt << " logdet: " << ll_logdet << " mpi: " << ll_mpi << std::endl;
        return -ll_sum; //invert for minimization
    }
};

//move to tests
int RealMatrixTest(void)
{
    RealMatrix Id = RealMatrix::Identity(3);
    RealMatrix B = RealMatrix::Zero(3, 3);
    RealMatrix CV(3, 1);
    CV(0, 0) = 1; CV(1, 0) = 2; CV(2, 0) = 3;
    RealMatrix RV(1, 3);
    RV(0, 0) = 1; RV(0, 1) = 2; RV(0, 2) = 3;


    //test RealMatrix::Diag()
    RealMatrix CV1(3, 1);
    CV1(0, 0) = 2; CV1(1, 0) = 4; CV1(2, 0) = 3;
    RealMatrix RM1(3, 3);
    RM1(0,0)=2;RM1(0,1)=0;RM1(0,2)=0;RM1(1,0)=0;RM1(1,1)=4;RM1(1,2)=0;RM1(2,0)=0;RM1(2,1)=0;RM1(2,2)=3;
    assert(RM1 == RealMatrix::Diag(CV1));

    assert(Id == Id*Id);
    assert(Id+Id == 2*Id);
    assert(Id-Id == B);
    assert(CV == Id*CV);
    assert(CV == RV.Transposed());
    assert(real(RV*CV) == real(14));

    assert(MatrixExponential(2*Id, 0) == Id);
    assert(MatrixExponential(2*Id, 1) == 3*Id);


    VectorMatrix IdV(2, 2);
    IdV(0, 0) = Vector(1, 2); IdV(0, 1) = Vector(0, 0);
    IdV(1, 0) = Vector(0, 0); IdV(1, 1) = Vector(3, 4);

    assert(IdV == IdV.Transposed());

    std::cout << Id.GetRows() << " " << Id.GetColumns() << std::endl;

    return 0;
}

int main(int argc, char** argv)
{
    CMDParameterEntropy args(argc, argv);

    real param;
    std::ofstream out;

    out.open(args.out_file);
    //out.precision(10);
    out << args.header() << std::endl;

    std::vector<unsigned int> times = MakeTimes(args.time_start, args.time_stop, args.time_step);

    Posdict positions;
    Angdict angles;

    //reads positions of all particles in data.txt into the vector data of real Vectors
    if (args.eval_model == Model::DynamicEntropy || args.eval_model == Model::DynamicOptimization || args.eval_model == Model::GeneralizedEntropy || args.eval_model == Model::GeneralizedEntropyJ)
    {
        read_data_consecutive(args.in_file, positions, angles, times);
    }
    else if (args.eval_model == Model::StaticEntropy)
    {
        read_data_single(args.in_file, positions, angles, times);
    }
    else
    {
        throw std::runtime_error("Invalid evaluation model.");
    }

    Pi pi = Pi(angles);
    if (pi.size() <= 0)
    {
        throw std::runtime_error("Data was not ready successfully.");
    }

    std::cout << "Initialization complete" << std::endl;

    if (args.eval_model == Model::DynamicEntropy)
    {
        std::vector<real> ni(args.sim_N);
        MATRIX n(args.sim_N, args.sim_N);
        n.Zero();
        out << "t param empnc JV(J) eta(T) c1s cs gs cts gts cint cpint gint chs ctint Likelihood" << std::endl;

        for (unsigned int t : times)
        {
            for (param = args.param_start; param < args.param_stop; param += args.param_step)
            {
                switch (args.eval_geometry)
                {
                case Geometry::topological:
                    calculate_n_topological(n, positions[t], args.sim_L, int(param));
                    assert(is_symmetric(n));
                    //calculate_n_topological_unsymmetrical(n, positions[t], args.sim_L, int(param));
                    //assert(is_number_of_neighbours_equal(n));
                    break;
                case Geometry::metric:
                    calculate_n_metric(n, positions[t], args.sim_L, param);
                    assert(is_symmetric(n));
                    break;
                case Geometry::voronoi:
                    calculate_n_voronoi(n, positions[t], args.sim_L);
                    assert(is_symmetric(n));
                    break;
                }

                create_ni(n, ni);

                real empnc = average(ni); //equivalent to get_average_number_of_neighbours(n); //empirical

                real c1s = C1s(t, pi);
                real cs = Cs(t, pi);
                real gs = Gs(t, pi);
                real cts = CTs(t, pi, empnc, ni);
                real gts = GTs(t, pi, empnc, ni);
                real cint = Cint(t, pi, empnc, n);
                real cpint = CPint(t, pi, empnc, n);
                real gint = Gint(t, pi, empnc, n);
                real chs = CHs(t, pi, empnc, ni);
                real ctint = CTint(t, pi, empnc, ni, n);

                real omega = Omega(args.sim_deltat, cint, gint);
                real temperature0 = T0(args.sim_deltat, cs, gs);
                real temperature0t = T0t(args.sim_deltat, cts, gts);

                real interaction = J(empnc, temperature0t, omega, cpint, chs, ctint);
                real temperature = T(args.sim_deltat, empnc, interaction, temperature0, omega, c1s, cs, gs, cts, gts, cpint, chs, ctint);
                real etaemp = eta_of_t(temperature);
                real JVemp = JV_of_J(args.sim_deltat, interaction, empnc);

                real likelihood = DynamicLogLikelihoodMaximized(args.sim_N, c1s, cs, gs, cts, gts, cint, cpint, gint, chs, ctint);

                //output
                std::cout << t << " " << param << " " << empnc<< " " << JVemp << " " << etaemp << " " << likelihood << std::endl;
                out << t << " " << param << " " << empnc << " " << JVemp << " " << etaemp << " " << c1s << " " << cs << " " << gs << " " << cts << " " << gts << " " << cint << " " << cpint << " " << gint << " " << chs << " " << ctint << " " << likelihood << std::endl;
            }
        }
    }
    else if (args.eval_model == Model::StaticEntropy)
    {
        std::vector<real> ni(args.sim_N);
        MATRIX n(args.sim_N, args.sim_N);
        n.Zero();
        out << "t param empnc jsigma2ratio cts cint ev_logsum ev0counter Likelihood" << std::endl;

        //real det;
        unsigned int r, c;
        unsigned int eval_num, ev0counter;
        real ev_logsum, eval_logsum, likelihood, jsigma2_ratio;
        for (unsigned int t : times)
        {
            for (param = args.param_start; param < args.param_stop; param += args.param_step)
            {
                switch (args.eval_geometry)
                {
                case Geometry::topological:
                    calculate_n_topological(n, positions[t], args.sim_L, int(param));
                    assert(is_symmetric(n));
                    //calculate_n_topological_unsymmetrical(n, positions[t], args.sim_L, int(param));
                    //assert(is_number_of_neighbours_equal(n));
                    break;
                case Geometry::metric:
                    calculate_n_metric(n, positions[t], args.sim_L, param);
                    assert(is_symmetric(n));
                    break;
                case Geometry::voronoi:
                    calculate_n_voronoi(n, positions[t], args.sim_L);
                    assert(is_symmetric(n));
                    break;
                }

                create_ni(n, ni);

                real empnc = average(ni); //equivalent to get_average_number_of_neighbours(n); //empirical

                real cts = CTs(t, pi, empnc, ni);
                real cint = Cint(t, pi, empnc, n);

                EigenMatrixReal lambda(args.sim_N, args.sim_N);

                for (r = 0; r < n.GetHeight(); ++r)
                {
                    for (c = 0; c < n.GetWidth(); ++c)
                    {
                        if (r == c)
                        {
                            //n_i delta_ij - n_ij
                            lambda(r, c) = ni[r] - real(*n.Get(r, c));
                        }
                        else
                        {
                            lambda(r, c) = - real(*n.Get(r, c));
                        }
                    }
                }

                //\lambda + \lambda^T
                //EigenMatrixReal lambdalambdaT = lambda + lambda.transpose();

                assert(is_weakly_diagonally_dominant_and_non_negativ_diagonal_entries(lambda));

                //eigenvalues of \lambda + \lambda^T
                Eigen::SelfAdjointEigenSolver<EigenMatrixReal> es;
                es.compute(lambda, /* computeEigenvectors = */ Eigen::EigenvaluesOnly);
                assert(es.info() == Eigen::ComputationInfo::Success);

                auto eval_vector = es.eigenvalues();

                //std::cout << "The eigenvalues of lambda are: " << ev_vector.transpose() << std::endl;

                eval_logsum = 0;
                ev_logsum = 0;
                ev0counter = 0;
                jsigma2_ratio = 1./(2.*empnc*(cts-cint));
                for (eval_num = 0; eval_num < eval_vector.size(); ++eval_num)
                {
                    assert(eval_vector(eval_num) > -0.0001);
                    if (eval_vector(eval_num) >= 0.0001)
                    {
                        //std::cout << ev_vector(evn) << " ";
                        // eval_logsum += std::log(jsigma2_ratio*eval_vector(eval_num)/PI);
                        ev_logsum += std::log(eval_vector(eval_num));
                        eval_logsum += std::log(jsigma2_ratio/PI);
                    }
                    else
                    {
                        ev0counter += 1;
                    }
                }
                //std::cout << lambda << std::endl;
                //std::cout << std::endl;

                eval_logsum += ev_logsum;

                likelihood = eval_logsum / 2. - args.sim_N / 2.;

                //output
                std::cout << t << " " << param << " " << empnc << " " << jsigma2_ratio << " " << ev0counter << " " << likelihood << std::endl;
                out << t << " " << param << " " << empnc << " " << jsigma2_ratio << " " << cts << " " << cint << " " << ev_logsum << " " << ev0counter << " " << likelihood << std::endl;
            }
        }
    }
    else if (args.eval_model == Model::DynamicOptimization)
    {
        dlib_vector starting_point(2); //3

        out << "time n_c empnc J sigma2 Likelihood" << std::endl;

        for (unsigned int t : times)
        {
            BufferedN n(args.sim_L, positions[t]);

            switch (args.eval_geometry)
            {
            case Geometry::voronoi:
                n.InitVoronoi();
                break;
            case Geometry::metric:
                n.InitMetric(args.param_start, args.param_stop, args.param_step);
                break;
            case Geometry::topological:
                n.InitTopological(args.param_start, args.param_stop, args.param_step);
                break;
            }

            n.Calculate();

            VectorMatrix pi0 = convert_to_vector_colvector(pi.Get(t));
            VectorMatrix pi1 = convert_to_vector_colvector(pi.Get(t+1));

            DynamicLikelihood l(args.sim_deltat, n, pi0, pi1);

            for (param = args.param_start; param < args.param_stop; param += args.param_step)
            {
                l.set_param(param);

                starting_point = real(0.75), real(0.2);

                //dlib::bfgs_search_strategy()
                //dlib::find_min_using_approximate_derivatives(
                //    dlib::cg_search_strategy(),
                //    dlib::objective_delta_stop_strategy(1e-7).be_verbose(), //14
                //    l,
                //    starting_point,
                //    -1
                //);

                dlib::matrix<real> lower_bound_constraint(2,1);
                lower_bound_constraint = -0.5, -0.5; // J, sigma2
                dlib::matrix<real> upper_bound_constraint(2,1);
                upper_bound_constraint = 4, 1; // J, sigma2

                try
                {
                    dlib::find_min_bobyqa(l,
                        starting_point,
                        6,    // number of interpolation points
                        lower_bound_constraint,  // lower bound constraint
                        upper_bound_constraint,   // upper bound constraint
                        0.1,    // initial trust region radius
                        1e-6,  // stopping trust region radius
                        1000    // max number of objective function evaluations
                    );
                }
                catch (const dlib::bobyqa_failure &f)
                {
                    std::cout << "Error for " << StringVector(starting_point) << ": " << f.what() << std::endl;
                    out << t << " " << param << " " << n.empnc(param) << " 0. 0. 0." << std::endl;
                    continue;
                }
                std::cout << t << " " << param << " " << n.empnc(param) << " " << StringVector(starting_point) << " " << -l(starting_point) << std::endl;
                out << t << " " << param << " " << n.empnc(param) << " " << StringVector(starting_point) << " " << -l(starting_point) << std::endl;
            }
        }
    }
    // optimizes L for each data point individually
    else if (args.eval_model == Model::GeneralizedEntropySingle)
    {
        out << "time n_c empnc J sigma2 Likelihood" << std::endl;

        for (unsigned int t : times)
        {
            BufferedLambda lambda(args.sim_L, positions[t]);

            switch (args.eval_geometry)
            {
                case Geometry::voronoi:
                    lambda.InitVoronoi();
                    break;
                case Geometry::metric:
                    lambda.InitMetric(args.param_start, args.param_stop, args.param_step);
                    break;
                case Geometry::topological:
                    lambda.InitTopological(args.param_start, args.param_stop, args.param_step);
                    break;
            }

            lambda.Calculate();

            VectorMatrix pi0 = convert_to_vector_colvector(pi.Get(t));
            VectorMatrix pi1 = convert_to_vector_colvector(pi.Get(t+1));

            GeneralLikelihoodSingle l(args.sim_deltat, lambda, pi0, pi1);

            dlib_vector starting_point(2); //3
            for (param = args.param_start; param < args.param_stop; param += args.param_step)
            {
                l.set_param(param);

                starting_point = real(0.75), real(0.2);

                //dlib::bfgs_search_strategy()
                //dlib::find_min_using_approximate_derivatives(
                //    dlib::cg_search_strategy(),
                //    dlib::objective_delta_stop_strategy(1e-7).be_verbose(), //14
                //    l,
                //    starting_point,
                //    -1
                //);

                dlib::matrix<real> lower_bound_constraint(2,1);
                lower_bound_constraint = -0.5, -0.5; // J, sigma2
                dlib::matrix<real> upper_bound_constraint(2,1);
                upper_bound_constraint = 4, 1; // J, sigma2

                try
                {
                    dlib::find_min_bobyqa(l,
                        starting_point,
                        6,    // number of interpolation points
                        lower_bound_constraint,  // lower bound constraint
                        upper_bound_constraint,   // upper bound constraint
                        0.1,    // initial trust region radius
                        1e-6,  // stopping trust region radius
                        1000    // max number of objective function evaluations
                    );
                }
                catch (const dlib::bobyqa_failure &f)
                {
                    std::cout << "Error for " << StringVector(starting_point) << ": " << f.what() << std::endl;
                    out << t << " " << param << " " << lambda.empnc(param) << " 0. 0. 0." << std::endl;
                    continue;
                }
                std::cout << t << " " << param << " " << lambda.empnc(param) << " " << StringVector(starting_point) << " " << -l(starting_point) << std::endl;
                out << t << " " << param << " " << lambda.empnc(param) << " " << StringVector(starting_point) << " " << -l(starting_point) << std::endl;
            }
        }
    }
    // optimizes L for sum of data points
    else if (args.eval_model == Model::GeneralizedEntropy)
    {

        out << "time n_c empnc J sigma2 Likelihood" << std::endl;

        GeneralLikelihoodSumWithJApprox l(args.sim_deltat, times, pi);

        for (param = args.param_start; param < args.param_stop; param += args.param_step)
        {
            l.Init(positions, param, args.eval_geometry, args.sim_L);
            l.CalculateUUI();
            real empnc = average(l.empncs);

            dlib_vector starting_point(2); //3

            starting_point = real(0.75), real(0.2);

            //dlib::bfgs_search_strategy()
            //dlib::find_min_using_approximate_derivatives(
            //    dlib::cg_search_strategy(),
            //    dlib::objective_delta_stop_strategy(1e-7).be_verbose(), //14
            //    l,
            //    starting_point,
            //    -1
            //);

            dlib::matrix<real> lower_bound_constraint(2,1);
            lower_bound_constraint = -0.5, -0.5; // J, sigma2
            dlib::matrix<real> upper_bound_constraint(2,1);
            upper_bound_constraint = 20, 8; // J, sigma2 4, 1

            try
            {
                dlib::find_min_bobyqa(l,
                    starting_point,
                    6,    // number of interpolation points
                    lower_bound_constraint,  // lower bound constraint
                    upper_bound_constraint,   // upper bound constraint
                    0.1,    // initial trust region radius
                    1e-6,  // stopping trust region radius
                    1000    // max number of objective function evaluations
                );
            }
            catch (const dlib::bobyqa_failure &f)
            {
                std::cout << "Error for " << StringVector(starting_point) << ": " << f.what() << std::endl;
                out << 0 << " " << param << " " << empnc << " 0. 0. 0." << std::endl;
                continue;
            }

            std::cout << 0 << " " << param << " " << empnc << " " << StringVector(starting_point) << " " << -l(starting_point) << std::endl;
            out << 0 << " " << param << " " << empnc << " " << StringVector(starting_point) << " " << -l(starting_point) << std::endl;
        }
    }
    // optimizes L for sum of data points, without J approximation
    else if (args.eval_model == Model::GeneralizedEntropyJ)
    {

        out << "time n_c empnc J sigma2 Likelihood" << std::endl;

        GeneralLikelihoodSumWithoutJApprox l(args.sim_deltat, times, pi);

        for (param = args.param_start; param < args.param_stop; param += args.param_step)
        {
            l.Init(positions, param, args.eval_geometry, args.sim_L);
            l.CalculateUUI();
            real empnc = average(l.empncs);

            dlib_vector starting_point(2); //3

            starting_point = real(0.75), real(0.2);

            //dlib::bfgs_search_strategy()
            //dlib::find_min_using_approximate_derivatives(
            //    dlib::cg_search_strategy(),
            //    dlib::objective_delta_stop_strategy(1e-7).be_verbose(), //14
            //    l,
            //    starting_point,
            //    -1
            //);

            dlib::matrix<real> lower_bound_constraint(2,1);
            lower_bound_constraint = -0.5, -0.5; // J, sigma2
            dlib::matrix<real> upper_bound_constraint(2,1);
            upper_bound_constraint = 20, 8; // J, sigma2 4, 1

            try
            {
                dlib::find_min_bobyqa(l,
                    starting_point,
                    6,    // number of interpolation points
                    lower_bound_constraint,  // lower bound constraint
                    upper_bound_constraint,   // upper bound constraint
                    0.1,    // initial trust region radius
                    1e-6,  // stopping trust region radius
                    1000    // max number of objective function evaluations
                );
            }
            catch (const dlib::bobyqa_failure &f)
            {
                std::cout << "Error for " << StringVector(starting_point) << ": " << f.what() << std::endl;
                out << 0 << " " << param << " " << empnc << " 0. 0. 0." << std::endl;
                continue;
            }

            std::cout << 0 << " " << param << " " << empnc << " " << StringVector(starting_point) << " " << -l(starting_point) << std::endl;
            out << 0 << " " << param << " " << empnc << " " << StringVector(starting_point) << " " << -l(starting_point) << std::endl;
        }
    }

    std::cout << "Evaluation complete" << std::endl;

    out.close();

    return 0;
}
