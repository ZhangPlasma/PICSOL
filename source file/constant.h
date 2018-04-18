#pragma once

/* physical constant */
#define EPS_0 8.85418782e-12  	// F/m, vacuum permittivity
#define MU_0 1.256637061e-6     // N/A^2, vacuum permeability 
#define K	1.38065e-23			// J/K, Boltzmann constant
#define ME 9.10938215e-31		// kg, electron mass
#define QE 1.602176565e-19		// C, elementary charge
#define AMU  1.660538921e-27	// kg, atomic mass unit

/* simulation flag */
#define GATHER_METHOD 1         // 0 for NGP, 1 for CIC
#define BOUND_CONDITION 0       // 0 for Dirichlet, 1 for Neumann, 2 for periodic

/* simulation constant */
#define SOR 1.6                 // parameter for SOR method, [0,2] 
#define LBOUND 0                // left boundary for x
#define RBOUND 1                // right boundary for x
#define MAX_LOOP 1000           // maximum cycle number
#define GRID_NUM 1000           // grid point number
#define MAX_PARTICLE 100000     // maximum particle number
#define MAX_ITERATION 100000    // maximum iteration for Gauss-Seidel method
#define MIN_ERROR  1e-4         // minimum error tolerance for iteration method    
#define ITERATION_CHECK 50      // iteration number to check convergence
#define RECORD_CHECK 100        // iteration number to 
#define DELTA_T 0.1             // minimum time step, normalized to inverse of plasma frequency 1/omega_{pe}
#define DELTA_X 1               // minimum space step, normalized to Debye length 

/* plasma parameters */
#define ETEMP 20                // electron temperature in eV
#define ITEMP 10                // ion temperature in eV
#define PLASMA_DEN   1.0e18     // average plasma density in m^-3
#define DEU 2.014102            // AU, deuterium mass

