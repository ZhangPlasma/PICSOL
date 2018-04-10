#include "Field.h"

Field::Field(int grid_num, double xbound1, double xbound2)
{
	if (grid_num <= 3)
	{
		cerr << " GRID NUMBER MUST BE A POSITIVE INTEGER GREATER THAN THREE!" << endl;
		system("pause");

	}
	this->grid_num = grid_num;
	this->xbound1 = xbound1;
	this->xbound2 = xbound2;
	this->inc = (xbound2 - xbound1) / (grid_num - 1);
	val = new double[grid_num];
	memset(val, 0, sizeof(double)* grid_num);
}

Field::~Field()
{
	delete[] val;
}

void Field::saveInfo(string directory)
{
	ofstream output_file(directory, ios::out);
	if (!output_file.is_open())
	{
		cerr << " FIELD OUTPUT FILE OPEN ERROR !" << endl;
		system("pause");
	}
	for (int i = 0; i < grid_num; i++)
	{
		output_file << val[i] << endl;
	}
}

void Field::clear()
{
	memset(val, 0, sizeof(double)* grid_num);
}

void Field::multiply(double constant)
{
	for (int i = 0; i < grid_num; i++)
	{
		val[i] *= constant;
	}
}

double Field::mean()
{
	double sum = 0;
	for (int i = 0; i < grid_num; i++)
	{
		sum += val[i];
	}
	return sum / grid_num;
}