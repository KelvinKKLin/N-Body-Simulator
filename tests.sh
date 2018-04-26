#!/bin/bash

rm results.txt
mpirun -n 2 ./x.project 1000 1000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 4 ./x.project 1000 1000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 8 ./x.project 1000 1000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 16 ./x.project 1000 1000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 32 ./x.project 1000 1000 0 9 1 0.25 2550 1960 ima >> results.txt

echo "\n\n\n\n" >> results.txt
mpirun -n 2 ./x.project 2000 2000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 4 ./x.project 2000 2000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 8 ./x.project 2000 2000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 16 ./x.project 2000 2000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 32 ./x.project 2000 2000 0 9 1 0.25 2550 1960 ima >> results.txt

echo "\n\n\n\n" >> results.txt
mpirun -n 2 ./x.project 4000 4000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 4 ./x.project 4000 4000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 8 ./x.project 4000 4000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 16 ./x.project 4000 4000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 32 ./x.project 4000 4000 0 9 1 0.25 2550 1960 ima >> results.txt

echo "\n\n\n\n" >> results.txt
mpirun -n 2 ./x.project 8000 1000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 4 ./x.project 8000 1000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 8 ./x.project 8000 1000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 16 ./x.project 8000 1000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 32 ./x.project 8000 1000 0 9 1 0.25 2550 1960 ima >> results.txt

echo "\n\n\n\n" >> results.txt
mpirun -n 2 ./x.project 16000 16000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 4 ./x.project 16000 16000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 8 ./x.project 16000 16000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 16 ./x.project 16000 16000 0 9 1 0.25 2550 1960 ima >> results.txt
mpirun -n 32 ./x.project 16000 16000 0 9 1 0.25 2550 1960 ima >> results.txt
