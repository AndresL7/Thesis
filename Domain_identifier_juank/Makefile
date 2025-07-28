CC = mpicc

#CFLAGS =  -c -O3 -I/home/$(USER)/local/include/ -I/usr/include/
#LFLAGS = -lm -L/home/$(USER)/local/lib -Wl,-R /home/$(USER)/local/lib

CFLAGS = -c -g -I/home/$(USER)/qhull/include/ -I/usr/include/ -I/usr/include/hdf5/serial/
LFLAGS = -lm -L/home/$(USER)/qhull/lib/ -Wl,"-R /home/$(USER)/qhull/lib/" -L/usr/lib/x86_64-linux-gnu/hdf5/serial

PROGRAM = Domain_identifier

$(PROGRAM):
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) $(OPENMP) -lqhull_r -lgsl -lgslcblas -lm -lhdf5 -o $@


clean:
	rm -rf $(PROGRAM)
	rm -rf *~
	rm -rf *#
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm -rf *.x
plot:
	gnuplot script.gpl
