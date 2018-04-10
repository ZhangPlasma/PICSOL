#include "algorithmn.h"

void gather(Species& particle, Field& Density, Field& Current)
{	
	/* constant parameter */
	static const double weight_den = PLASMA_DEN / (MAX_PARTICLE / GRID_NUM);
	double xbound = Density.xbound1;
	double inc = Density.inc;
	
	/* field value reset to zero */
	Density.clear();
	Current.clear();

	/* field gathering process */
	if (GATHER_METHOD == 0)  /* NGP(Nearest Grid Point) method implement */
	{ 
		/* update density and current field */
		for (int i = 0; i < particle.part.size(); i++)
		{   
			/* nearest grid method */
			int pos = round((particle.part[i].x1 - xbound) / inc);
			Density.val[pos] ++;
			Current.val[pos] += particle.part[i].v1;
		}
	}
	else if (GATHER_METHOD == 1)  /* CIC(Cloud In Cell) method implement */
	{
		/* update density and current field */
		for (int i = 0; i < particle.part.size(); i++)
		{
			double ratio = (particle.part[i].x1 - xbound) / inc;
			int pos = floor(ratio);
			double prop = ratio - pos;
            
			/* linear interploration */
			Density.val[pos] += 1 - prop;	
			Density.val[pos + 1] += prop;
			Current.val[pos] += particle.part[i].v1 * (1 - prop);
			Current.val[pos + 1] += particle.part[i].v1 * prop;
		}
	}
	else
	{
		cerr << " GATHER METHOD ERROR !" << endl;
		system("pause");
	}

    /* normaliztion for physical variables */
	for (int i = 0; i < Density.grid_num; i++)
	{
		Density.val[i] *= weight_den;
		Current.val[i] *= particle.charge * weight_den;
	}
}

void rhoComp(Field& electDen, Field& ionDen, Field& rho) 
{
	for (int i = 0; i < rho.grid_num; i++)
	{
		rho.val[i] = QE * (ionDen.val[i] - electDen.val[i]);
	}
}

int poissonSOR(Field& phi, Field& rho)
{
	/* pre-computed coefficients */
	int grid_num = phi.grid_num;
	double dx2 = phi.inc * phi.inc;
	phi.val[0] = phi.val[grid_num - 1] = 0;   /* Dirichlet boundaries */

	/* solve potential using iterative method */
	for (int solver_iter = 0; solver_iter < MAX_ITERATION; solver_iter++)
	{
		for (int i = 1; i < grid_num - 1; i++)
		{   
			double phi_prime = 0.5*(phi.val[i - 1] + phi.val[i + 1] + dx2 * rho.val[i]);
			phi.val[i] = phi.val[i] + SOR * (phi_prime - phi.val[i]);
		}

		/* check for convergence */
		if (solver_iter % ITERATION_CHECK == 0)
		{
			double sum = 0;
			for (int i = 1; i < grid_num - 1; i++)
			{   
				double R = - rho.val[i] - (phi.val[i - 1] - 2 * phi.val[i] + phi.val[i + 1]) / dx2;
				sum += R * R;
			}
			
			if ((sqrt(sum) / grid_num) < MIN_ERROR) 
			{
				return solver_iter; 
			}
		}
	}
	cerr << "POISSON SOLVER CONVERGENCE ERROR !" << endl;
	system("pause");
}

void poissonMatrix(Field& phi, Field& rho)
{
	/* pre-computed coefficients */
	int grid_num = phi.grid_num;
	double dx2 = phi.inc * phi.inc;
	static const double a = 1;
	static const double b = -2;
	static const double c = 1;
	double* c_prime = new double[grid_num - 3]; /* Dirichlet boundaries */
    
    /* Dirichlet boundaries */
	phi.val[0] = 0;
	phi.val[grid_num - 1] = 0;
	
	/* use phi.val to store d_prime[grid_num -2] */
	for (int i = 1; i < grid_num - 1; i++)
	{
		phi.val[i] = - rho.val[i] * dx2 / EPS_0;
	}
    
	/* iteration intializer */
	c_prime[0] = c / b;
	phi.val[1] /= b;

	/* modify the coefficients forward */
	for (int i = 1; i < grid_num - 2; i++)
	{
		double temp = b - c_prime[i - 1] * a;
		c_prime[i] = c / temp; 
		phi.val[i + 1] = (phi.val[i + 1] - phi.val[i] * a) / temp;
	}

	/* back substitute from d_prime to phi */
	for (int i = grid_num - 3; i > 0; i--)
	{
		phi.val[i] = phi.val[i] - c_prime[i - 1] * phi.val[i + 1];
	}	
}

void electFieldComp(Field& phi, Field& electField)
{  
	/* central difference */
	for (int i = 1; i < electField.grid_num - 1; i++)
	{
		electField.val[i] = (phi.val[i - 1] - phi.val[i + 1]) / (2 * electField.inc);
	}
    
	/* forward difference at boundaries */
	electField.val[0] = (phi.val[0] - phi.val[1]) / electField.inc;
	electField.val[electField.grid_num - 1] = (phi[phi.grid_num - 2] - phi[phi.grid_num - 1]) / electField.inc;
}

void pusher(Species& particle, Field& electField, double delta_t)
{
	/* particle pusher process must be in pair with field gather process */
	double xbound1 = electField.xbound1;
	double xbound2 = electField.xbound2;
	double inc = electField.inc;
	double sc = particle.specific_charge;

	if (GATHER_METHOD == 0)  /* NGP(Nearest Grid Point) method implement */
	{
		for (int i = 0; i < particle.part.size(); i++)
		{
			/* particle's electric field is the same as the nearset grid point */
			int pos = round((particle.part[i].x1 - xbound1) / inc);

			/* advance particle */
			particle.part[i].x1 += particle.part[i].v1 * delta_t;
			particle.part[i].v1 += sc * electField.val[pos] * delta_t;

			if (particle.part[i].x1 < xbound1 || particle.part[i].x1 > xbound2)
			{
				/* particle collides with boundary */
				particle.part.erase(particle.part.begin() + i);
			}
		}
	}
	else if (GATHER_METHOD == 1)  /* CIC(Cloud In Cell) method implement */
	{	
		for (int i = 0; i < particle.part.size(); i++)
		{
			/* particle's electric field is the linear interploration of near grid point */
			double ratio = (particle.part[i].x1 - xbound1) / inc;
			int pos = floor(ratio);
			double prop = ratio - pos;

			/* advance particle */
			particle.part[i].x1 += particle.part[i].v1 * delta_t;
			particle.part[i].v1 += sc * delta_t * ((1 - prop) * electField.val[pos] + prop * electField.val[pos + 1]);

			if (particle.part[i].x1 < xbound1 || particle.part[i].x1 > xbound2)
			{
				/* particle collides with boundary */
				particle.part.erase(particle.part.begin() + i);
			}
		}	
	}
	else
	{
		cerr << " PUSHER METHOD ERROR !" << endl;
		system("pause");
	}
}

void leapFrogRewind(Species& particle, Field& electField, double delta_t)
{
	double xbound1 = electField.xbound1;
	double xbound2 = electField.xbound2;
	double inc = electField.inc;
	double sc = particle.specific_charge;

	if (GATHER_METHOD == 0)  /* NGP(Nearest Grid Point) method implement */
	{
		for (int i = 0; i < particle.part.size(); i++)
		{
			/* particle's electric field is the same as the nearset grid point */
			int pos = round((particle.part[i].x1 - xbound1) / inc);

			/* velocity rewinding */
			particle.part[i].v1 -= 0.5 * sc * electField.val[pos] * delta_t;
		}
	}
	else if (GATHER_METHOD == 1)  /* CIC(Cloud In Cell) method implement */
	{
		for (int i = 0; i < particle.part.size(); i++)
		{
			/* particle's electric field is the linear interploration of near grid point */
			double ratio = (particle.part[i].x1 - xbound1) / inc;
			int pos = floor(ratio);
			double prop = ratio - pos;

			/* velocity rewinding */
			particle.part[i].v1 -= 0.5 * sc * delta_t * ((1 - prop) * electField.val[pos] + prop * electField.val[pos + 1]);
		}
	}
	else
	{
		cerr << " REWINDING METHOD ERROR !" << endl;
		system("pause");
	}
}







