
rm libtest.so
rm *.o *.lo
mpicc -c -acc -fPIC -Wfatal-errors -I. create_xgrid_gpu.c create_xgrid_utils_gpu.c general_utils_gpu.c  grid_utils.c tree_utils.c
mpicc -acc -shared create_xgrid_gpu.o create_xgrid_utils_gpu.o general_utils_gpu.o grid_utils.o tree_utils.o -o libtest.so
