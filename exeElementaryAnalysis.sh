g++ NewtonianSpace.hpp
g++ -c -o NewtonianSpace.o NewtonianSpace.cpp
g++ -c -o ElementaryAnalysis.o ElementaryAnalysis.cpp
g++ -o exeElementaryAnalysis NewtonianSpace2.o ElementaryAnalysis.o
./exeElementaryAnalysis
