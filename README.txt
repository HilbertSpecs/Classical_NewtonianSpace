				MAC OSX Terminal

Instructions for Compiling and Running ElementaryAnalysis2.cpp source code file.
 		
Class Specification file: NewtonianSpace.hpp
Class Implementation file: NewtonianSpace2.cpp
Source Code: ElementaryAnalysis2.cpp

Object Files created: g++ -c -o NewtonianSpace2.o NewtonianSpace2.cpp
			g++ -c -o ElementaryAnalysis.o ElementaryAnalysis.cpp


Link the Two Object files created: NewtonianSpace2.o , ElementaryAnalysis.o
			g++ -o exe2ElementaryAnalysis NewtonianSpace2.o ElementaryAnalysis.o

Executable Created: exe2ElementaryAnalysis

Run Program:	./exe2ElementaryAnalysis

Program will Print to following Output Files:	GeneratedAcceleration.txt
						CalculatedForce.txt
						PositionTrajectory.txt
						VelocityTrajectory.txt
						WorkFunction.txt
						PowerFunction.txt

CopyRight: Evan James Rabeaux
5/12/2017
This Software was written in full by Evan James Rabeaux.
Please do not redistribute, copy, or use this software for any other use than as a sample of my coding work.

