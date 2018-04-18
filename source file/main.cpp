#include "stdafx.h"
#include "constant.h"
#include "Field.h"
#include "Particle.h"
#include "Species.h"
#include "algorithm.h"

/**************************************************************************************************** 
 PICSOL is a 1d-1v elctrostatic(2d-3V electromagnetic in the future) PIC( Particle In Cell) 
 code developed by Huaxiang Zhang since 2018 Mar.28th, which is primarily designed to 
 simulate SOL(Scrape-off Layer) plasma in tokamaks.This code is open-source to anyone who is 
 intended to use it. The latest version can be downloaded from my github homepage 
 https://github.com/ZhangPlasma/PICSOL/
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

	/* create data storage file based on current system time */
	time_t t = time(NULL);
	struct tm current_time;
	localtime_s(&current_time, &t);
	string directory = string("./data/") + to_string(current_time.tm_year) + string("-") 
		                                 + to_string(current_time.tm_mon) + string("-")
		                                 + to_string(current_time.tm_mday) + string(" ")
		                                 + to_string(current_time.tm_hour) + string(":")
		                                 + to_string(current_time.tm_min) + string("/");
/***********************************************************************************************************************************/
/***********************************************************************************************************************************/
    
	/* MAIN LOOP */
	start = clock();
	for (int main_iter = 0; main_iter < MAX_LOOP; main_iter++)
	{	
		/* particle pusher */
		pusher(electron, electField, delta_t);
		pusher(ion, electField, delta_t);

		/* field update */
		gather(electron, electDen, electCurrent);
		gather(ion, ionDen, ionCurrent);
		rhoComp(electDen, ionDen, rho);
	    poissonMatrix(phi, rho);                      
		electFieldComp(phi, electField);

		/* save data */
		if (main_iter % RECORD_CHECK == 0)
		{
			/* save field data */
			electCurrent.saveInfo(directory + to_string(main_iter) + string("electCurrent.dat"));
			electDen.saveInfo(directory + to_string(main_iter) + string("electDen.dat"));
			ionCurrent.saveInfo(directory + to_string(main_iter) + string("ionCurrent.dat"));
			ionDen.saveInfo(directory + to_string(main_iter) + string("ionDen.dat"));
			phi.saveInfo(directory + to_string(main_iter) + string("phi.dat"));
			rho.saveInfo(directory + to_string(main_iter) + string("rho.dat"));
			electField.saveInfo(directory + to_string(main_iter) + string("electField.dat"));

			/* save particle data */
			electron.saveParticle(directory + string("electron.dat"));
			ion.saveParticle(directory + string("ion.dat"));
		}
	}
	end = clock();
	cout << "Main Loop is complete, for " << MAX_LOOP << " steps, time consuming: " << 1000 * (end - start) / (double)(CLOCKS_PER_SEC) << "ms" << endl;

	system("pause");
	return 0;
}
