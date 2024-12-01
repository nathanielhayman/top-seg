#!bin/bash
gcc -shared -DDEBUG -I ./CBLAS/include remove.c -L./CBLAS/lib -lcblas -L./CBLAS/Open_BLAS/lib -lopenblas -o libremove.dll -fPIC
