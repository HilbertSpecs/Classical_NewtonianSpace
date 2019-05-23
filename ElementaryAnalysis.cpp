#include <iostream>
#include <cmath>
#include <iomanip>
#include "NewtonianSpace.hpp"
using namespace std;

/*
 Parameter in object is the size of the dynamically allocated arrays for
 radial position, particle velocity, acceleration, and Force. It must match the
 1st parameter in generatedAcceleration().
*/


int main()
{


	NewtonianSpace object(100);

	object.setMass(40.5);
	object.setIntPosition(0 ,0,0);
	object.setIntVelocity(0,0,0);
	object.setIntAcceleration(0,0,0);
	object.getMass();
	object.getIntPosition();
	object.getIntVelocity();
	object.getIntAcceleration();
	object.generateAcceleration(100,1,20);
	object.calculateForce(100);	
	object.computeVelocity(100, .005);
	object.computePosition(100, .005);
	object.velocityTrajectory(100);
	object.positionTrajectory(100);
	object.computeWork(100);
	object.computePower(100, .005);
	object.calculateWork(100);
	object.calculatePower(100);

return 0;

}






/*

CopyRight: Evan James Rabeaux
5/12/2017
This Software was written in full by Evan James Rabeaux.
Please do not redistribute, copy, or use this software for any other use than as a viewing sample of my coding work.

*/