rm -f *.txt
make
mpirun -n $1 dist/app