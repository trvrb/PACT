# Compiling for Unix: make
# Compiling for Windows: make CROSS=i386-mingw32-

CC=$(CROSS)g++
LD=$(CROSS)ld
AR=$(CROSS)ar

pact: main.o node.o coaltree.o series.o io.o param.o rng.o
	$(CC) -O3 -o pact main.o node.o coaltree.o series.o io.o param.o rng.o
main.o: main.cpp node.h coaltree.h series.h io.h param.h rng.h
	$(CC) -O3 -c main.cpp 
node.o: node.cpp node.h 
	$(CC) -O3 -c node.cpp 
coaltree.o: coaltree.cpp coaltree.h 
	$(CC) -O3 -c coaltree.cpp 
series.o: series.cpp series.h 
	$(CC) -O3 -c series.cpp 	
io.o: io.cpp io.h 
	$(CC) -O3 -c io.cpp 	
param.o: param.cpp param.h 
	$(CC) -O3 -c param.cpp 
rng.o: rng.cpp rng.h 
	$(CC) -O3 -c rng.cpp 	
clean: 
	rm *.o pact
