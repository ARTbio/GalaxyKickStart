/*
 * main.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: rocha
 */

#include "./lib/out/studentttests.h"
#include <iostream>
#include <fstream>

namespace Support
{
   static const unsigned int MAX_READS = 1024*256;
   static const unsigned int INVALID_IDENTIFICATION = -1;
   static const char END_OF_INFORMATION[2] = ";";
   static char SEPARATOR[2] = ",";
   static const char END_OF_LINE[2] = "\0";
   static const char COMMENTARY[2] = "#";
   static const double CONFIDENCE_FACTOR = 0.95;
}

static std::vector < std::vector < double > > * readCSVfile(const char * fileName)
{
	std::ifstream file(fileName, std::ios::in);
	char in[Support::MAX_READS];
	char value[Support::MAX_READS];
	std::vector < std::vector < double > > * table = new std::vector < std::vector < double > > (0);
	std::vector < double > sample(0);
	unsigned int subPosition, bytesRead, position;

	if (not file)
	{
		std::cerr << fileName << " not found!" << std::endl;
		return NULL; //exit(1);
	}

	do
		file.getline(in, Support::MAX_READS);
	while (in[0] == Support::COMMENTARY[0] and not file.eof());
	if (file.eof()) return NULL;
	bytesRead = file.gcount();
	while(bytesRead >= 2)
	{
		table->push_back(sample);
		position = 0;
		do
		{
			subPosition = 0;
			do
				value[subPosition++] = in[position++];
			while (position < bytesRead and in[position] != Support::SEPARATOR[0] and in[position] != Support::END_OF_INFORMATION[0]);
			value[subPosition] = Support::END_OF_LINE[0];
			table->back().push_back(atof(value));
			position++;
			subPosition = 0;
		}
		while (position < bytesRead and in[position -1] != Support::END_OF_INFORMATION[0]);

		do file.getline(in, Support::MAX_READS);
		while (in[0] == Support::COMMENTARY[0] and not file.eof());
		if (file.eof()) break;
		bytesRead = file.gcount();
	}
	file.close();

	return table;
}

int main(int numberOfArguments, char ** arguments)
{
	ap::real_1d_array x, y;
	std::vector < std::vector < double > > * xTable, * yTable;
    double bothtails = 0, lefttail = 0, righttail = 0;
    Support::SEPARATOR[2] = arguments[3][0];
    unsigned int xColumn = atoi(arguments[4]);
    unsigned int yColumn;

    if (numberOfArguments > 5) yColumn = atoi(arguments[5]);
    else yColumn = xColumn;

    xTable = readCSVfile(arguments[1]);
    yTable = readCSVfile(arguments[2]);

    x.setlength(xTable->size());
	y.setlength(yTable->size());

	for (unsigned int xElement=0; xElement < xTable->size(); xElement++)
		x(xElement) = xTable->at(xElement).at(xColumn);
	for (unsigned int yElement=0; yElement < yTable->size(); yElement++)
		y(yElement) = yTable->at(yElement).at(yColumn);

	studentttest2(x, xTable->size(), y, yTable->size(), bothtails, lefttail, righttail);

	std::cout << bothtails << "," << lefttail << "," << righttail << std::endl;

	delete xTable;
	delete yTable;

	return 0;
}
