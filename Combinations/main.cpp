/*
 * main.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: fabiorjvieira
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char ** args)
{
	unsigned int n = atoi(args[1]);
	bool in[n], oneZero, zeroOne;
	char elements[n];

	for (unsigned int i = 0; i < n; i++)
	{
		in[i] = false;
		elements[i] = '1' + i;
	}

	for (unsigned j = 0; j < pow(2,n); j++)
	{
		oneZero = false;
		zeroOne = false;
		for (unsigned int i = 1; i < n-1; i++)
		{
			if (in[i-1] == true and in[i] == false)
				oneZero = true;
			if (in[i] == false and in[i+1] == true and oneZero)
				zeroOne = true;
		}

		for (int i = n-1; i >= 0; i--)
			std::cout << in[i];
		std::cout << " - ";

		if (oneZero and zeroOne)
		{
			for (unsigned int i = 0; i < n; i++)
				if (in[i] == true)
					std::cout << elements[i] << ",";
		}
		std::cout << std::endl;

		for (unsigned int i = 0; i < n; i++)
		{
			if (in[i] == false)
			{
				in[i] = true;
				break;
			}
			else in[i] = false;
		}
	}

	return 0;
}
