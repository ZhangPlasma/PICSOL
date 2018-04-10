#include "stdafx.h"
#include "constant.h"
#include "Field.h"
#include "Particle.h"
#include "Species.h"
#include "algorithmn.h"

/**************************************************************************************************** 
 PICSOL is a 1d-1v elctrostatic(2d-3V electromagnetic in the future) PIC( Particle In Cell) 
 code developed by Huaxiang Zhang since 2018 Mar.28th, which is primarily designed to 
 simulate SOL(Scrape-off Layer) plasma in tokamaks.This code is open-source to anyone who is 
 intended to use it. The latest version can be downloaded from my github homepage 

 Version Number: 1.01
******************************************************************************************************/

int main()
{   

	/* RNG seeding */
	srand((unsigned int)time(NULL));

	/* time and space normalization */
	double delta_t = DELTA_T;
	double delta_x = DELTA_X;

	/* particle initialization */
	clock_t start = clock();
	Species electron(ME, QE);
	Species ion(AMU*DEU, QE);
	electron.initialMaxwell(MAX_PARTICLE, ETEMP, LBOUND, RBOUND);
	ion.initialMaxwell(MAX_PARTICLE, ITEMP, LBOUND, RBOUND);
	clock_t end = clock();
	cout << "Particle initialization is complete, time consuming: " << 1000 * (end - start) / (double)(CLOCKS_PER_SEC) << "ms" << endl;

	/* field declare */
	start = clock();
	Field electField(GRID_NUM, LBOUND, RBOUND);
	Field phi(GRID_NUM, LBOUND, RBOUND);
	Field electDen(GRID_NUM, LBOUND, RBOUND);
	Field ionDen(GRID_NUM, LBOUND, RBOUND);
	Field rho(GRID_NUM, LBOUND, RBOUND);
	Field electCurrent(GRID_NUM, LBOUND, RBOUND);
	Field ionCurrent(GRID_NUM, LBOUND, RBOUND);

	/* field initialization */
	gather(electron, electDen, electCurrent);
	gather(ion, ionDen, ionCurrent);
	rhoComp(electDen, ionDen, rho);
	poissonMatrix(phi, rho);  
	electFieldComp(phi, electField);

	/* leapfrog rewinding */
	leapFrogRewind(electron, electField, delta_t);
	leapFrogRewind(ion, electField, delta_t);
	end = clock();
	cout << "Field initialization is complete, time consuming: " << 1000 * (end - start) / (double)(CLOCKS_PER_SEC) << "ms" << endl;

/***********************************************************************************************************************************/
/***********************************************************************************************************************************/
    
	/* MAIN LOOP */
	start = clock();
	for (int main_iter = 0; main_iter < MAX_LOOP; main_iter++)
	{	
		/* particle pusher */
		pusher(electron, electField, delta_t);
		pusher(ion, electField, delta_t);

		/* field updating */
		gather(electron, electDen, electCurrent);
		gather(ion, ionDen, ionCurrent);
		rhoComp(electDen, ionDen, rho);
	    poissonMatrix(phi, rho);                      
		electFieldComp(phi, electField);

		/* record data */
		if (main_iter % RECORD_CHECK == 0)
		{
			/* create data storage file based on current system time */
			time_t t = time(NULL);
			char currentTime[64];
			strftime(currentTime, sizeof(currentTime), "%Y-%m-%d %H-%M-%S", localtime(&t));
			string directory("./data/");
			directory += string(currentTime);

			/* record field data */
			electCurrent.saveInfo(directory + string("electCurrent.dat"));
			electDen.saveInfo(directory + string("electDen.dat"));
			ionCurrent.saveInfo(directory + string("ionCurrent.dat"));
			ionDen.saveInfo(directory + string("ionDen.dat"));
			phi.saveInfo(directory + string("phi.dat"));
			rho.saveInfo(directory + string("rho.dat"));
			electField.saveInfo(directory + string("electField.dat"));

			/* record particle data */
			electron.saveParticle(directory + string("electron.dat"));
			ion.saveParticle(directory + string("ion.dat"));
		}
        
	}
	end = clock();
	cout << "Main Loop is complete, for " << MAX_LOOP << " steps, time consuming: " << 1000 * (end - start) / (double)(CLOCKS_PER_SEC) << "ms" << endl;

	system("pause");
	return 0;
}
