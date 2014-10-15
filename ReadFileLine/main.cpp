#include <stdlib.h>
#include <fstream>
#include <iostream>

void readFileLine(const char * fileName, unsigned int lineNumber)
{
	unsigned int maxRead = 1024000;
	char in[maxRead];
	std::ifstream file(fileName, std::ios::in);

	lineNumber--;
	for (unsigned int line = 0; line < lineNumber; line++)
		file.ignore(maxRead,'\n');
	file.getline(in,maxRead);
	std::cout << in << std::endl;
}


int main(int nargs, char ** args)
{
	readFileLine(args[1],atoi(args[2]));
	return 0;
}
