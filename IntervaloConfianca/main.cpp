/*
 * main.cpp
 *
 *  Created on: Jun 5, 2010
 *      Author: rocha
 */

#include "main.hpp"

double confidence(double confidenceFactor, double std, unsigned int size)
{
   return (invstudenttdistribution(size, confidenceFactor) * std)/std::sqrt((double)size);
}

double standardDerivation(std::vector < double > values, double average)
{
   if (values.size() == 1) return 0;
   double std = 0;

   for (unsigned int element = 0; element < values.size(); element++)
      std += std::pow(average - values[element],2);
   std /= (values.size() -1);

   return std::sqrt(std);
}

static double ** readCSVfile(const char * fileName)
{
	std::ifstream file(fileName, std::ios::in);
	char in[Support::MAX_READS];
	char value[Support::MAX_READS];
	unsigned int subPosition, columnSize = 0;
	unsigned int tableLineQuantity = 0, tableColumnQuantity = 0;
	double ** table;

	unsigned int bytesRead, position;

	if (not file)
	{
		std::cerr << fileName << " not found!" << std::endl;
		return NULL; //exit(1);
	}

	/*comments*///computing maximum column size and line quantity
	file.getline(in, Support::MAX_READS);
	bytesRead = file.gcount();
	while(bytesRead > 3)
	{
		if (in[0] != Support::COMMENTARY[0])
		{
			position = 0;
			columnSize = 0;
			do
			{
				subPosition = 0;
				do
					value[subPosition++] = in[position++];
				while (position < bytesRead and in[position] != Support::SEPARATOR[0] and in[position] != Support::END_OF_INFORMATION[0]);
				position++;
				subPosition = 0;
				columnSize++;
			}
			while (position < bytesRead and in[position -1] != Support::END_OF_INFORMATION[0]);
			tableLineQuantity++;
			if (tableColumnQuantity < columnSize) tableColumnQuantity = columnSize;
		}
		if (not file.eof())
		{
			file.getline(in, Support::MAX_READS);
			bytesRead = file.gcount();
		}
		else break;
	}
	file.close();

	/*comments*///Allocation of routes
	table = new double * [tableLineQuantity +1];
	for (unsigned int line = 0; line <= tableLineQuantity; line++)
	{
		table[line] = new double [tableColumnQuantity +1];
		for (unsigned int column = 0; column <= tableColumnQuantity; column++)
			table[line][column] = Support::INVALID_IDENTIFICATION;
	}

	/*comments*///load routes
	file.open(fileName, std::ios::in);
	for (unsigned int line = 0; line < tableLineQuantity; line++)
	{
		do file.getline(in, Support::MAX_READS);
		while (in[0] == Support::COMMENTARY[0]);
		bytesRead = file.gcount();

		position = 0;
		columnSize=0;
		do
		{
			subPosition = 0;
			do
				value[subPosition++] = in[position++];
			while (position < bytesRead and in[position] != Support::SEPARATOR[0] and in[position] != Support::END_OF_INFORMATION[0]);
			value[subPosition] = Support::END_OF_LINE[0];
			table[line][columnSize] = atof(value);
			columnSize++;
			position++;
			subPosition = 0;
		}
		while (position < bytesRead and in[position -1] != Support::END_OF_INFORMATION[0]);
	}
	file.close();

	return table;
}

void intervaloConfianca(const char * inFileName, const char * outFileName, unsigned int column)
{
	double ** table;
	std::string fileName;
	std::vector < double > values;
	double standardDerivationValue, confidenceValue, average, total=0;
	unsigned int line;

    std::cout << "Reading:" << inFileName << " ";
    std::cout.flush();
    table = readCSVfile(inFileName);

    if (table != NULL)
    {
    	for (line = 0; table[line][0] != Support::INVALID_IDENTIFICATION; line++)
    	{
    		total += table[line][column];
    		values.push_back(table[line][column]);
    		delete [] table[line];
    	}

    	std::cout << "Read:" << line << std::endl;
    	delete [] table[line];
    	delete [] table;
    }

	std::cout << "Creating:" << outFileName << std::endl;
	std::ofstream file(outFileName, std::ios::out);
	if (not file)
	{
		std::cerr << fileName << " not created!" << std::endl;
		exit(1);
	}
	file << Support::COMMENTARY << "Performance:media,devio,confianca,limite sup, limite inf" << std::endl;

	if (values.size() != 0)
	{
		average = total / values.size();
		standardDerivationValue = standardDerivation(values, average);
		confidenceValue = confidence(Support::CONFIDENCE_FACTOR, standardDerivationValue, values.size());
		file << average << Support::SEPARATOR << standardDerivationValue << Support::SEPARATOR << confidenceValue << Support::SEPARATOR << (average + confidenceValue) << Support::SEPARATOR << (average - confidenceValue) << Support::END_OF_INFORMATION << std::endl;
	}

	file.close();
}

int main(int numberOfArguments, char ** arguments)
{
	int result = 0;

	intervaloConfianca(arguments[1],arguments[2],atoi(arguments[3]));

	return result;
}
