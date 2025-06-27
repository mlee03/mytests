CC = mpicc

clib.so: c.o
	$(CC) -shared $^ -o clib.so

c.o: c.c
	$(CC) -fPIC -c $^
