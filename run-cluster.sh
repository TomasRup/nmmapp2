rm -f *.txt
mpicc -std=c99 -L/usr/lib/x86_64-linux-gnu -lm -o dist/app app/main.c
sbatch mpi-sbatch