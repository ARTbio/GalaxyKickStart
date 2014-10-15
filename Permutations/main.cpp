/*
 * main.cpp
 *
 *  Created on: Apr 4, 2014
 *      Author: fabiorjvieira
 */

#include <vector>
#include <iostream>
#include <algorithm>

void runPermutations(std::vector < unsigned long int > permutations, std::vector < unsigned long int > options)
{
	unsigned int currentOption = 0;
	while (options.size() > 0)
	{
		permutations.push_back(options.at(currentOption));
		options.erase(options.begin() + currentOption);
		currentOption++;
		runPermutations(permutations, options);
		if (options.size() == 0)
		{
			for (unsigned int permutationIndex = 0; permutationIndex < permutations.size(); permutationIndex++)
				std::cout << " " << permutations.at(permutationIndex);
			std::cout << std::endl;
		}
	}
}

int main(int argn, char ** argc)
{
	std::vector < unsigned long int > permutations, options;

	for (unsigned int optionIndex = 0; optionIndex < 4; optionIndex++)
		options.push_back(optionIndex);

	 std::sort (myints,myints+3);

	  std::cout << "The 3! possible permutations with 3 elements:\n";
	  do {
	    std::cout << myints[0] << ' ' << myints[1] << ' ' << myints[2] << '\n';
	  } while ( std::next_permutation(myints,myints+3) );

	  std::cout << "After loop: " << myints[0] << ' ' << myints[1] << ' ' << myints[2] << '\n';


	runPermutations(permutations, options);

	return 0;
}
