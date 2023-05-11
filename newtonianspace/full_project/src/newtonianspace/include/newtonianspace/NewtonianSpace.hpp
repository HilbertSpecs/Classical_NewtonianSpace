#pragma once

#ifndef NEWTONIANSPACE
#define NEWTONIANSPACE
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cmath>
using namespace std;


const double DEFAULT = 0.0;

class NewtonianSpace
{

		//When creating classes I will use all the variables in the equation as my private variables as a minimum.

	private:

		int n, min, max, size;
		double dt;
		double *x, *y, *z, *rmag;
		double *vx, *vy, *vz, *vmag;
		double *ax, *ay, *az, *amag;
		double *Fx, *Fy, *Fz,*Fmag, m;
		double *Wx, *Wy, *Wz, *Wmag, WxTotal, WyTotal, WzTotal, WmagTotal;
		double *Px, *Py, *Pz, *Pmag, PxTotal, PyTotal, PzTotal, PmagTotal;

		void create_xArray(int size)
		{
			n = size;
		 	x = new double[size];
		for( int i = 0; i < size; i++)
			x[i] = DEFAULT;
		}
		
		void create_yArray(int size)
		{
			n = size;
			y = new double[size];
		for(int i= 0; i < size; i++)
			y[i] = DEFAULT;
		}
		
		void create_zArray(int size)
		{
			n = size;
			z = new double[size];
		for(int i = 0; i < size;i++)
			z[i] = DEFAULT;
		}
		
		void create_rmagArray(int size)
		{
			n = size;
			rmag = new double[size];
		for(int i=0; i < size; i++)
			rmag[i] = DEFAULT;
		}

		void create_vxArray(int size)
		{
			n= size;
			vx = new double[size];
		for(int i= 0; i < size; i++)
			vx[i] = DEFAULT;
		}

		void create_vyArray(int size)
		{
			n= size;
			vy = new double[size];
		for(int i =0; i < size; i++)
			vy[i] = DEFAULT;
		}

		void create_vzArray(int size)
		{
			n = size;
			vz = new double[size];
		for(int i =0; i < size; i++)
			vz[i] = DEFAULT;
		}

		void create_vmagArray(int size)
		{
			n = size;
			vmag = new double[size];
		for(int i = 0; i < size; i++)
			vmag[i] = DEFAULT;
		}

		void create_axArray(int size)
		{
			n = size;
			ax = new double[size];
		for(int i = 0; i < size; i++)
			ax[i] = DEFAULT;
		}
		
		void create_ayArray(int size)
		{
			n = size;
			ay = new double[size];
		for(int i = 0; i < size; i++)
			ay[i] = DEFAULT;
		}

		void create_azArray(int size)
		{
			n= size;
			az = new double[size];
		for(int i = 0; i < size; i++)
			az[i] = DEFAULT;
		}

		void create_amagArray(int size)
		{
			n = size;
			amag = new double[size];
		for( int i = 0; i < size; i++)
			amag[i] = DEFAULT;
		}

		void create_FxArray(int size)
		{
			n = size;
			Fx = new double[size];
		for( int i =0; i < size; i++)
			Fx[i] = DEFAULT;
		}

		void create_FyArray(int size)
		{
			n = size;
			Fy = new double[size];
		for(int i = 0; i < size; i++)
			Fy[i]= DEFAULT;
		}

		void create_FzArray(int size)
		{
			n = size;
			Fz = new double[size];
		for(int i = 0 ; i < size; i++)
			Fz[i] = DEFAULT;
		}
		
		void create_FmagArray(int size)
		{
			n = size;	
			Fmag = new double[size];
		for(int i =0; i < size; i++)
			Fmag[i] = DEFAULT;
		}		


		void create_WxArray(int size)
		{
			n = size;
			Wx = new double[size];
		for(int i= 0; i < size; i++)
			Wx[i] = DEFAULT;
		}

		void create_WyArray(int size)
		{
			n = size;
			Wy = new double[size];
		for(int i =0; i < size; i++)
			Wy[i] = DEFAULT;
		}
		
		void create_WzArray(int size)
		{
			n = size;
			Wz = new double[size];
		for(int i = 0; i < size; i++)
			Wz[i] = DEFAULT;
		}

		void create_WmagArray(int size)
		{
			n = size;
			Wmag = new double[size];
		for(int i = 0; i <size ; i++)
			Wmag[i] = DEFAULT;
		}

		void create_PxArray(int size)
		{
			n = size;
			Px = new double[size];
		for(int i= 0; i < size; i++)
			Px[i] = DEFAULT;
		}

		void create_PyArray(int size)
		{
			n = size;
			Py = new double[size];
		for(int i= 0; i < size; i++)
			Py[i] = DEFAULT;
		}

		void create_PzArray(int size)
		{
			n= size;
			Pz = new double[size];
		for(int i= 0; i < size; i++)
			Pz[i] = DEFAULT;
		}

		void create_PmagArray(int size)
		{
			n = size;
			Pmag = new double[size];
		for(int i=0; i< size; i++)
			Pmag[i] = DEFAULT;
		}





		//When creating classes I will use all functions I want that can parametrize these private variables.
		// Constructors, Assigning or Setting Parameters, Getting these Values, Printing these values, Destructors.

	public:
        NewtonianSpace();
		NewtonianSpace(int);
		~NewtonianSpace();
		virtual void setMass(double mass);
		void setIntPosition(double, double, double);
		void setIntVelocity(double, double, double);
		void setIntAcceleration(double, double, double);
		virtual double getMass() const;
		double getIntPosition() const;
		double getIntVelocity() const;
		double getIntAcceleration() const;
		void generateAcceleration(int,int,int);		
		void calculateForce(int);
		void computeVelocity(int, double);
		void computePosition(int, double);
		void velocityTrajectory(int);
		void positionTrajectory(int);
		void computeWork(int);
		void computePower(int, double);	
		double calculateWork(int);
		double calculatePower(int);
};

#endif


/*
*****************************************************************************************************************************************
Written & Distributed by: RABOMETRICS
*****************************************************************************************************************************************
*/
