pact: main.o node.o coaltree.o series.o tree.hh.gch
	g++ -o pact main.o node.o coaltree.o series.o
main.o: main.cpp node.h coaltree.h 
	g++ -c main.cpp 
node.o: node.cpp node.h 
	g++ -c node.cpp 
coaltree.o: coaltree.cpp coaltree.h 
	g++ -c coaltree.cpp 
series.o: series.cpp series.h 
	g++ -c series.cpp 	
tree.hh.gch: tree.hh
	g++ -c tree.hh
clean: 
	rm *.o pact *.gch