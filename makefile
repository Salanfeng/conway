all:
	gcc conway.c -std=c11 -o conway1 -fopenmp -Ofast -funroll-loops -march=native -mfma -mavx2 -m3dnow
