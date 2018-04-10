#pragma once
#include "stdafx.h"
#include "Constant.h"
#include "Field.h"
#include "Species.h"

/*
gather - Specialized function gathers field information on the grid from the particle distribution and external condition
*/
void gather(Species& particle, Field& Density, Field& Current);

/*
rhoComp - compute charge density 
*/
void rhoComp(Field& electDen, Field& ionDen, Field& rho);

/*
poissonSOR - 1D poisson solver using SOR(succsive overrelaxation) iterative method 
*/
int poissonSOR(Field& phi, Field& rho);

/*
poissonMatrix - 1D poisson solver using Liewellyn Thomas's algorithmn for solving tridiagonal matrix
*/
void poissonMatrix(Field& phi, Field& rho);

/*
electFieldComp - compute electric field using electrostatic potential 
*/
void electFieldComp(Field& phi, Field& electField);

/*
pusher - Specialized function move the particle based on the field quantities on the grid point and external condition
*/
void pusher(Species& particle, Field& electField, double delta_t);

/*
leapFrogRewind - Specialized function to rewind particle velocity to form a staggered leap frog format
*/
void leapFrogRewind(Species& particle, Field& electField, double delta_t);
