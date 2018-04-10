#include "Species.h" 

void Species::initialMaxwell(int num, double temp, double xbound1, double xbound2) 
{   
	/* uniformlly distributed in x space, maxwell distributed in v space */
	for (int i = 0; i < num; i++)
	{
		double x = xbound1 + (xbound2 - xbound1)*rand() * (1.0 / RAND_MAX);
		double vt = sqrt(temp * QE / mass);
		double v = generateMaxwell(0, vt);
		part.push_back(Particle<double>(x, v));
	}
}

void Species::saveParticle(string directory)
{
	ofstream output_file(directory, ios::out);
	if (!output_file.is_open())
	{
		cerr << " PARTICLE OUTPUT FILE OPEN ERROR !" << endl;
		system("pause");
	}
	for (int i = 0; i < (int)part.size(); i++)
	{
		output_file << part[i].x1 << part[i].v1 << endl;
	}
}

double Species::kinetic()
{
	double kinetic = 0;
	for (int i = 0; i < part.size(); i++)
	{
		kinetic += 1 / 2 * mass * sqrt(part[i].v1);
	}
	return kinetic;
}

void Species::addParticle(double x, double v)
{
	part.push_back(Particle<double>(x, v));
}

void Species::remParticle(int index)
{
	part.erase(part.begin() + index);
}