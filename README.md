# Vicsek Maximum likelihood estimation

tested with gcc version 4.9.2 (x86_64-win32-seh-rev2, Built by MinGW-W64 project)   
compiled with "-std=c++11", "-D_USE_MATH_DEFINES", "-DNDEBUG"  
tested with "-fexpensive-optimizations", "-O2", "-march=native"   
on Windows 7 e.g. `g++.exe -std=c++11  -fexpensive-optimizations -O2 -march=native -DNDEBUG -IC:\voro++-0.4.6\src -IC:\eigen-eigen-10219c95fe65 -IC:\dlib-18.16  -c entropy.cpp -o entropy.o`  
`g++.exe -std=c++11  -fexpensive-optimizations -O2 -march=native -DNDEBUG  -IC:\voro++-0.4.6\src -IC:\eigen-eigen-10219c95fe65 -IC:\dlib-18.16  -c MyMatrix.cpp -o MyMatrix.o`  
`g++.exe  -o EntropyVicsek.exe constants.o entropy.o MyMatrix.o tests.o C:\libVoro++.a`  


additional dependencies:
- [voro++](http://math.lbl.gov/voro++/)
- [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [dlib](http://dlib.net/)


call executable as `EntropyVicsek.exe data data_out staticentropy N L dt geometry geometry_start geometry_end geometry_step tstart tstop tstep`

where  
data (data): input file "data.txt"  
data_out (data_out): output file "data_out.txt"  
mode (dynamicentropy): staticentropy: [static entropy](http://www.pnas.org/content/109/13/4786.full), dynamicentropy: [dynamic entropy](http://arxiv.org/pdf/1310.3810v1.pdf), generalizedentropy: numeric optimisation, generalizedentropyJ: numeric optimisation without n_i \approx <n_i> restriction  
N (256): number of particles  
L (16.): length of square simulation box  
dt (0.01): time difference between consecutive snapshots  
geometry (metric): metric, topological, voronoi  
geometry_start (0.1): metric: lowest radius to calculate likelihood for, topological: lowest number of neighbours to calculate likelihood for, leave this parameter out for voronoi  
geometry_end (1.5): metric: highest radius, topological: highest number of neighbours, leave this parameter out for voronoi  
geometry_step (0.1): metric: radius step, topological: number of neighbours step, leave this parameter out for voronoi  
tstart (1): start from this line in input file data (beware additional header lines and factor 2 because of positions and angles)  
tstop (200): end in this line in input file data  
tstep (2): 1: snapshots a,b,c,d,... are evaluated a->b, b->c, ..., 2: evaluated a->b, c->d  
