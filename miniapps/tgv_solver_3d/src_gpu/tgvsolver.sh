
#gpu调试
# 编译 CPU 模块，生成调试信息
nvfortran -g -O0 -c commtype.F90 commvar.F90 constdef.F90 commfunc.F90  singleton.F90 filter.F90 derivative.F90

# 编译 GPU 模块（保留 device 调试信息）
nvfortran -gpu=debug,cc89 -g -lineinfo -O0 -c commarray_gpu.cuf commvar_gpu.cuf commfunc_gpu.cuf fludyna_gpu.cuf bc_gpu.cuf filter_gpu.cuf derivative_gpu.cuf solver_gpu.cuf

# 链接生成可执行文件（链接时同样保留调试信息）
nvfortran -gpu=debug,cc89 -g -O0 -I. tgvsolver.cuf \
    commtype.o commvar.o constdef.o commfunc.o  singleton.o filter.o derivative.o\
    commarray_gpu.o commvar_gpu.o commfunc_gpu.o fludyna_gpu.o bc_gpu.o filter_gpu.o derivative_gpu.o solver_gpu.o\
    -o tgvsolver
#调试运行
cuda-gdb ./tgvsolver
compute-sanitizer ./tgvsolver

#end gpu调试

#nsys调试
    nvfortran -O3 -g -c commtype.F90 commvar.F90 constdef.F90 commfunc.F90  singleton.F90 filter.F90 derivative.F90

    nvfortran -gpu=cc89 -O3 -g -lineinfo -c commarray_gpu.cuf commvar_gpu.cuf commfunc_gpu.cuf fludyna_gpu.cuf bc_gpu.cuf filter_gpu.cuf derivative_gpu.cuf solver_gpu.cuf

    nvfortran -gpu=cc89 -O3 -g -Mvect -I. tgvsolver.cuf \
    commtype.o commvar.o constdef.o commfunc.o  singleton.o filter.o derivative.o\
    commarray_gpu.o commvar_gpu.o commfunc_gpu.o fludyna_gpu.o bc_gpu.o filter_gpu.o derivative_gpu.o solver_gpu.o\
    -o tgvsolver

    nsys profile --force-overwrite true -o tgvsolver_nsys_report ./tgvsolver
    
    ncu --set full \
    --import-source yes \
    --source-folders . \
    --call-stack \
    --kernel-name-base demangled \
    -o tgvsolver_ncu_report \
    ./tgvsolver

    sudo /opt/nvidia/hpc_sdk/Linux_x86_64/26.1/compilers/bin/ncu --set full \
    --import-source yes \
    --source-folders . \
    --call-stack \
    -f \
    --kernel-name-base demangled \
    -o tgvsolver_ncu_report \
    ./tgvsolver

#end nsys调试

#正常运行
nvfortran -O3 -Kieee -c commtype.F90 commvar.F90 constdef.F90 commfunc.F90  singleton.F90 filter.F90 derivative.F90

nvfortran -O3 -gpu=cc89 -Kieee -c commarray_gpu.cuf commvar_gpu.cuf commfunc_gpu.cuf fludyna_gpu.cuf bc_gpu.cuf filter_gpu.cuf derivative_gpu.cuf solver_gpu.cuf

nvfortran -O3 -gpu=cc89 -Kieee -I. tgvsolver.cuf \
    commtype.o commvar.o constdef.o commfunc.o  singleton.o filter.o derivative.o\
    commarray_gpu.o commvar_gpu.o commfunc_gpu.o fludyna_gpu.o bc_gpu.o filter_gpu.o derivative_gpu.o solver_gpu.o\
    -o tgvsolver

./tgvsolver
    
#end正常运行



nvfortran -cuda -lineinfo -O0 -traceback -c commtype.F90 commvar.F90 constdef.F90 commfunc.F90  singleton.F90 filter.F90 derivative.F90

nvfortran -cuda -gpu=cc89 -lineinfo -O0 -traceback -c commarray_gpu.cuf commvar_gpu.cuf commfunc_gpu.cuf fludyna_gpu.cuf bc_gpu.cuf filter_gpu.cuf derivative_gpu.cuf solver_gpu.cuf

nvfortran -cuda -gpu=cc89 -O0 -traceback -I. tgvsolver.cuf \
    commtype.o commvar.o constdef.o commfunc.o  singleton.o filter.o derivative.o\
    commarray_gpu.o commvar_gpu.o commfunc_gpu.o fludyna_gpu.o bc_gpu.o filter_gpu.o derivative_gpu.o solver_gpu.o\
    -o tgvsolver

compute-sanitizer --tool memcheck ./tgvsolver



nvfortran -O3 -Kieee -I. tgvsolver.F90 \
    commtype.o commvar.o constdef.o commfunc.o  singleton.o filter.o derivative.o\
    -o tgvsolver