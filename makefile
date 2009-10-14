CC=$(CROSS)g++
LD=$(CROSS)ld
AR=$(CROSS)ar

pact: main.o node.o coaltree.o series.o io.o param.o
	$(CC) -o pact main.o node.o coaltree.o series.o io.o param.o
main.o: main.cpp node.h coaltree.h series.h io.h param.h
	$(CC) -c main.cpp 
node.o: node.cpp node.h 
	$(CC) -c node.cpp 
coaltree.o: coaltree.cpp coaltree.h 
	$(CC) -c coaltree.cpp 
series.o: series.cpp series.h 
	$(CC) -c series.cpp 	
io.o: io.cpp io.h 
	$(CC) -c io.cpp 	
param.o: param.cpp param.h 
	$(CC) -c param.cpp 
clean: 
	rm *.o pact
