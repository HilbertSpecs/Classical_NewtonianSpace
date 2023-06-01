//Need to go in to a few places and compute the rmag, vmag, amag, and Fmag vectors and print them to their according
// .txt files

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "include/newtonianspace/NewtonianSpace.hpp"
using namespace std;

NewtonianSpace::NewtonianSpace()
{
    size = 100;
    
    min=0;
    max=0;
    dt=0;
    m=0;
    WxTotal=0;
    WyTotal=0;
    WzTotal=0;
    WmagTotal=0;
    PxTotal=0;
    PyTotal=0;
    PzTotal=0;
    PmagTotal=0;

    create_xArray(size);
    create_yArray(size);
    create_zArray(size);
    create_rmagArray(size);
    create_vxArray(size);
    create_vyArray(size);
    create_vzArray(size);
    create_vmagArray(size);
    create_axArray(size);
    create_ayArray(size);
    create_azArray(size);
    create_amagArray(size);
    create_FxArray(size);
    create_FyArray(size);
    create_FzArray(size);
    create_FmagArray(size);
    create_WxArray(size);
    create_WyArray(size);
    create_WzArray(size);
    create_WmagArray(size);
    create_PxArray(size);
    create_PyArray(size);
    create_PzArray(size);
    create_PmagArray(size);

}

NewtonianSpace::NewtonianSpace(int len)
{
    min=0;
    max=0;
    dt=0;
    m=0;
    WxTotal=0;
    WyTotal=0;
    WzTotal=0;
    WmagTotal=0;
    PxTotal=0;
    PyTotal=0;
    PzTotal=0;
    PmagTotal=0;
    
    size = len;
    create_xArray(size);
    create_yArray(size);
    create_zArray(size);
    create_rmagArray(size);
    create_vxArray(size);
    create_vyArray(size);
    create_vzArray(size);
    create_vmagArray(size);
    create_axArray(size);
    create_ayArray(size);
    create_azArray(size);
    create_amagArray(size);
    create_FxArray(size);
    create_FyArray(size);
    create_FzArray(size);
    create_FmagArray(size);
    create_WxArray(size);
    create_WyArray(size);
    create_WzArray(size);
    create_WmagArray(size);
    create_PxArray(size);
    create_PyArray(size);
    create_PzArray(size);
    create_PmagArray(size);

}

NewtonianSpace::~NewtonianSpace()
{
	delete [] x;
	delete [] y;
	delete [] z;
	delete [] rmag;
	delete [] vx;
	delete [] vy;
	delete [] vz;
	delete [] vmag;
	delete [] ax;
	delete [] ay;
	delete [] az;
	delete [] amag;
	delete [] Fx;
	delete [] Fy;
	delete [] Fz;
	delete [] Fmag;
	delete [] Wx;
	delete [] Wy;
	delete [] Wz;
	delete [] Wmag;
	delete [] Px;
	delete [] Py;
	delete [] Pz;
	delete [] Pmag;
}

void NewtonianSpace::setMass(double mass)
{
	m = mass;
}

void NewtonianSpace::setIntPosition(double rx, double ry, double rz)
{
	x[0]= rx;
	y[0]= ry;
	z[0]= rz;	
	rmag[0] = sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));
}

void NewtonianSpace::setIntVelocity(double vox, double voy, double voz)
{
	vx[0] = vox;
	vy[0] = voy;
	vz[0] = voz;
	vmag[0] = sqrt(pow(vx[0],2)+pow(vy[0],2)+pow(vz[0],2));
}

void NewtonianSpace::setIntAcceleration(double aox, double aoy, double aoz)
{
	ax[0] = aox;
	ay[0] = aoy;
	az[0] = aoz;
	amag[0] = sqrt(pow(ax[0],2)+pow(ay[0],2)+pow(az[0],2));
}



double NewtonianSpace::getMass() const
{

	//cout << "The mass of the object in the Space is " << m << endl;

    return m;
}

double NewtonianSpace::getIntPosition() const
{
	//cout << "The initial position of the object relative to the origin is " << x[0] << ", " << y[0] << ", " << z[0] << ", " << rmag[0] << endl;

    /*return x[0];
    return y[0];
    return z[0];*/
    return rmag[0];

}

double NewtonianSpace::getIntVelocity() const
{
	//cout << "The initial velocity of the object is " << vx[0] << ", " << vy[0] << ", " << vz[0] << ", " << vmag[0] << endl;

    /*return vx[0];
    return vy[0];
    return vz[0];*/
    return vmag[0];

}

double NewtonianSpace::getIntAcceleration() const
{
	//cout << "The initial acceleration of the object is " << ax[0] << ", " << ay[0] << ", " << az[0] << ", " << amag[0] << endl;

    /*return ax[0];
    return ay[0];
    return az[0];*/
    return amag[0];

}

void NewtonianSpace::generateAcceleration(int n,int min,int max)
{
	fstream dataFile;
	int N, minimum, maximum;
	N = n;
	minimum = min;
	maximum = max;
 	
	dataFile.open("GeneratedAccel.txt", ios::out);
    
    dataFile << fixed << setprecision(4);
	dataFile << "The Acceleration values generated for " << N << " values are: " << endl;
	dataFile << "index\t\tax\t\tay\t\taz\t\tamag" << endl;
	dataFile << "0" << "\t\t" << ax[0] << "\t\t" << ay[0] << "\t\t" << az[0] << "\t\t" << amag[0] << endl;
 

	for(int i = 1; i < N; i++)
	{
		ax[i] = (rand()%maximum + minimum)*.37;
		ay[i] = (rand()%maximum + minimum)*.28;
		az[i] = (rand()%maximum + minimum)*.19;
		amag[i] = sqrt(pow(ax[i],2)+pow(ay[i],2)+pow(az[i],2));

	dataFile << i << "\t\t" << ax[i] << "\t\t" << ay[i] << "\t\t" << az[i] << "\t\t" << amag[i] << endl;

	}

    dataFile.close();

}

void NewtonianSpace::calculateForce(int size)
{	

	fstream dataFile;	
	int i,N;
	N = size;	
	Fx[0] = m * ax[0];
	Fy[0] = m * ay[0];
	Fz[0] = m * az[0];
	Fmag[0] = m * amag[0];

	dataFile.open("CalculatedForce.txt", ios::out);
	
    dataFile << fixed << setprecision(3);
	dataFile << "The computed Force for the " << N << " Generated Acceleration values are: " << endl;	
	dataFile << "index\t\tFx\t\tFy\t\tFz\t\tFmag" << endl;
	//cout << "0" << "\t\t" << Fx[0] << "\t\t" << Fy[0] << "\t\t" << Fz[0] << "\t\t" << Fmag[0] << endl;
	dataFile << "0" << "\t\t" << Fx[0] << "\t\t" << Fy[0] << "\t\t" << Fz[0] << "\t\t" << Fmag[0] << endl;


	for(int i= 1; i < N;i++)
		{
			Fx[i] = m*ax[i];
			Fy[i] = m*ay[i];
			Fz[i] = m*az[i];
			Fmag[i] = m*amag[i];	

		dataFile << i << "\t\t" << Fx[i] << "\t\t" << Fy[i] << "\t\t" <<  Fz[i] << "\t\t" << Fmag[i] << endl;


		}

    dataFile.close();

}

void NewtonianSpace::computeVelocity(int size, double timestep)
{

	int i, N;
	double dt;
	N = size;
	dt = timestep;
	for(int i = 1 ; i < N; i++)
	{	
		vx[i] += ax[i-1]*dt + vx[i-1];
		vy[i] += ay[i-1]*dt + vy[i-1];
		vz[i] += az[i-1]*dt + vz[i-1];
		vmag[i] = sqrt(pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2));

	}


}

void NewtonianSpace::computePosition(int size, double timestep)
{
	int i, N;
	double dt;
	N = size;
	dt = timestep;
	for(int i = 1; i < N; i++)
	{
		x[i] += vx[i-1]*dt + x[i-1];
		y[i] += vy[i-1]*dt + y[i-1];
		z[i] += vz[i-1]*dt + z[i-1];
		rmag[i] = sqrt(pow(x[i],2) + pow(y[i],2) + pow(z[i],2));

	}


}		

void NewtonianSpace::velocityTrajectory(int size)
{
	fstream dataFile;
	int N;
	N = size;
	
	dataFile.open("VelocityTrajectory.txt", ios::out);
	
    dataFile << fixed << setprecision(4);
	dataFile << "The Velocity Values for the object are: " << endl;
	dataFile << "index\t\tvx\t\tvy\t\tvz\t\tvmag" << endl;

	for(int i =0; i < N; i++)
	{
		dataFile << i << "\t\t" << vx[i] << "\t\t" << vy[i] << "\t\t" << vz[i] << "\t\t" << vmag[i] << endl;
	}

    dataFile.close();

}

void NewtonianSpace::positionTrajectory(int size)
{
	fstream dataFile;
	int N;
	N = size;
	
	dataFile.open("PositionTrajectory.txt", ios::out);

    dataFile << fixed << setprecision(4);
	dataFile << " The Position Values for the object are: " << endl;
	dataFile << "index\t\tx\t\ty\t\tz\t\trmag" << endl; 
	
	for(int i = 0; i < N; i++)
	{
		dataFile << i << "\t\t" << x[i] << "\t\t" << y[i] << "\t\t" << z[i] << "\t\t" << rmag[i] << endl;
	}

    dataFile.close();

}


void NewtonianSpace::computeWork(int size)
{
	fstream dataFile;
	int i, N;
	double *dx, *dy, *dz;
	N = size;
	
	dataFile.open("WorkFunction.txt", ios:: out);

	dx = new double[N];
	dy = new double[N];
	dz = new double[N];
	
	Wx[0] = 0;
	Wy[0] = 0;
	Wz[0] = 0;
    
    dataFile << fixed << setprecision(4);
    dataFile << "The Values of the Work as a Function of Time are: " << endl;
    dataFile << "index\t\tWx\t\tWy\t\tWz\t\tWmag" << endl;

	for(int i = 1 ; i < N; i++)
	{
		dx[i] = x[i] - x[i-1];
		dy[i] = y[i] - y[i-1];
		dz[i] = z[i] - z[i-1];
		
		Wx[i] = Fx[i]*dx[i];
		Wy[i] = Fy[i]*dy[i];
		Wz[i] = Fz[i]*dz[i];

		Wmag[i]= sqrt(pow(Wx[i],2)+pow(Wy[i],2)+pow(Wz[i],2));

		dataFile << i << "\t\t" << Wx[i] << "\t\t" << Wy[i] << "\t\t" << Wz[i] << "\t\t" << Wmag[i] << endl;

	}


    delete [] dx;
    delete [] dy;
    delete [] dz;


    dataFile.close();


}


void NewtonianSpace::computePower(int size, double timestep)
{
	fstream dataFile;
	int i, N;
	double dt;
	N = size;
	dt = timestep;	

    Px[0] = Wx[0]/dt;
    Py[0] = Wy[0]/dt;
    Pz[0] = Wz[0]/dt;
    Pmag[0] = sqrt(pow(Px[0],2)+pow(Py[0],2)+pow(Pz[0],2));
    
	dataFile.open("PowerFunction.txt", ios::out);	
    dataFile << fixed << setprecision(2);
    dataFile << "The Values for the Power as a function of Time are: " << endl;
    dataFile << "index\t\tPx\t\tPy\t\tPz\t\tPmag" << endl;
    dataFile << "0" << "\t\t" << Px[0] << "\t\t" << Py[0] << "\t\t" << Pz[0] << "\t\t" << Pmag[0] << endl;
    
    
	for( int i = 1; i < N; i++)
	{
		Px[i] = Wx[i]/dt;
		Py[i] = Wy[i]/dt;
		Pz[i] = Wz[i]/dt;
		Pmag[i] = sqrt(pow(Px[i],2)+pow(Py[i],2)+pow(Pz[i],2));

        dataFile << i << "\t\t" << Px[i] << "\t\t" << Py[i] << "\t\t" << Pz[i] << "\t\t" << Pmag[i] << endl;

	}
	

    dataFile.close();

}

double NewtonianSpace::calculateWork(int size)
{
	int i, N;
	N = size;
	
	for(int i = 0; i < N; i++)
	{
		WxTotal += Wx[i];
		WyTotal += Wy[i];
		WzTotal += Wz[i];
		WmagTotal += Wmag[i];
	}

	//cout << "WxTotal\t\tWyTotal\t\tWzTotal\t\tWmagTotal" << endl;
	//cout << WxTotal << "\t\t" << WyTotal << "\t\t" << WzTotal << "\t\t" << WmagTotal << endl;
	
    /*return WxTotal;
    return WyTotal;
    return WzTotal;*/
    return WmagTotal;

}

double NewtonianSpace::calculatePower(int size)
{
	int i, N;
	N = size;

	for(int i =0 ; i < N; i++)
	{
		PxTotal += Px[i];
		PyTotal += Py[i];
		PzTotal += Pz[i];
		PmagTotal += Pmag[i];
	}

	//cout << "PxTotal\t\tPyTotal\t\tPzTotal\t\tPmagTotal" << endl;
	//cout << PxTotal << "\t\t" << PyTotal << "\t\t" << PzTotal << "\t\t" << PmagTotal << endl;

    /*return PxTotal;
    return PyTotal;
    return PzTotal;*/
    return PmagTotal;

}

/*
*****************************************************************************************************************************************
Written & Distributed by: RABOMETRICS
*****************************************************************************************************************************************
*/



