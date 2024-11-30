#!bin/bash
gcc -shared -I ./CBLAS/include remove.c -L./CBLAS/lib -lcblas -L./CBLAS/Open_BLAS/lib -lopenblas -o libremove.dll -fPIC
