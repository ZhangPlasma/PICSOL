#pragma once
#include "stdafx.h"

/*
generateMaxwell - generate maxwell distridution random number N(mu,sigma) based on Box-Muller Method
*/
double generateMaxwell(double mu, double sigma)
{
	/* static constant to ensure double number boundary */
	static const double epsilon = numeric_limits<double>::min();
	static const double two_pi = 2.0*3.14159265358979323846;

	/* Box-Muller Method */
	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);

	double z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	return z0 * sigma + mu;
}
