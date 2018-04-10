#pragma once
#include "stdafx.h"
/*
Field - data structure to store field value
*/

class Field
{
public:

	/* variable list */
    int grid_num;        // number of grid point
	double* val;         // Field value storage
	double xbound1;      // left bound of grid region
	double xbound2;      // right bound of grid region 
	double inc;          // displasment between each grid point

	/* funtion list */
	Field(int grid_num, double xbound1, double xbound2);     // uniformlly spaced grid
	void saveInfo(string directory);                         // save gird value in the given directory
	void multiply(double constant);                          // grid value multiply by a constant
	double mean();                                           // average of field value  
	void clear();                                            // reset grid value to zero
	virtual ~Field();                                        // overwrite distructor function to avoid memory leak
};

