all: sssp_serial sssp_parallel

sssp_serial: src/util.h src/util.c src/sssp_serial.c
	gcc src/util.c src/sssp_serial.c -o sssp_serial
	
sssp_parallel: src/util.h src/util.c src/sssp_parallel.c
	mpicc src/util.c src/sssp_parallel.c -o sssp_parallel

clean:
	rm sssp_serial sssp_parallel out.txt out2.txt
