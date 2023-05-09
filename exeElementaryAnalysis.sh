g++ NewtonianSpace.hpp
g++ -c -o NewtonianSpace2.o NewtonianSpace2.cpp
g++ -c -o ElementaryAnalysis.o ElementaryAnalysis.cpp
g++ -o exeElementaryAnalysis NewtonianSpace2.o ElementaryAnalysis.o
./exeElementaryAnalysis
