CC = mpicc

#CFLAGS = -c -g -Wall -I/home/$(USER)/local/include/ -I/usr/include/
#LFLAGS = -lm -Wall -L/home/$(USER)/local/lib 

CFLAGS = $(OPTIONS) $(OPENMP) -c -g -I/home/$(USER)/local/include/ -I/usr/include/ -I/usr/include/hdf5/serial/
LFLAGS = -lm -L/home/$(USER)/local/lib -Wl,"-R /home/$(USER)/local/lib" -L/usr/lib/x86_64-linux-gnu/hdf5/serial

PROGRAM = Delaunay5

Delaunay5:
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o  $(LFLAGS) -lm -o $@.x -lqhull_r -lgsl -lgslcblas -lhdf5

clean:
	rm -rf $(PROGRAM)
	rm -rf *~
	rm -rf *.out
	rm -rf *.x
	rm -rf *#
	rm -rf *.o
	rm -rf *.a	
	rm -rf *.so
	rm -rf *.btr
