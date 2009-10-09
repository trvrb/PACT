pact: main.o node.o coaltree.o series.o io.o param.o
	g++ -o pact main.o node.o coaltree.o series.o io.o param.o
main.o: main.cpp node.h coaltree.h series.h io.h param.h
	g++ -c main.cpp 
node.o: node.cpp node.h 
	g++ -c node.cpp 
coaltree.o: coaltree.cpp coaltree.h 
	g++ -c coaltree.cpp 
series.o: series.cpp series.h 
	g++ -c series.cpp 	
io.o: io.cpp io.h 
	g++ -c io.cpp 	
param.o: param.cpp param.h 
	g++ -c param.cpp 
clean: 
	rm *.o pact