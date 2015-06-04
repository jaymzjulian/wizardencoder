g++ -DUSE_OPENMP -O6 -Wall -Iresid -Iinih encoder.cc resid/*.cc inih/*.c -o test -lfftw3 -fopenmp  fastsid/*.c fastsid/*.cc
#g++ -O6 -Wall -Iresid -Iinih encoder.cc resid/*.cc inih/*.c -o test -lfftw3 fastsid/*.c fastsid/*.cc
#g++ -O6 -Wall -Iresid -Iinih encoder.cc resid/*.cc inih/*.c -o test -lfftw3 
#g++ -g -Wall -Iresid -Iinih encoder.cc resid/*.cc inih/*.c -o test -lfftw3  fastsid/*.c fastsid/*.cc
