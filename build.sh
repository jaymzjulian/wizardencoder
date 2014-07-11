g++ -DUSE_OPENMP -O6 -Wall -Iresid -Iinih encoder.cc resid/*.cc inih/*.c -o test -lfftw3 -fopenmp 
#g++ -O6 -Wall -Iresid -Iinih encoder.cc resid/*.cc inih/*.c -o test -lfftw3 
#g++ -O0  -g -Wall -Iresid -Iinih encoder.cc resid/*.cc inih/*.c -o test -lfftw3 -fopenmp 
