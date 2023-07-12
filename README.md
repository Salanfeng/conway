gcc更改版本至10.2.0

编译(或直接用makefile):	
```
gcc conway.c -std=c11 -o conway1 -fopenmp -Ofast -funroll-loops -march=native -mfma -mavx2 -m3dnow
```
环境变量：
```
export OMP_PROC_BIND=true
```
运行 
```
sh run.sh
```