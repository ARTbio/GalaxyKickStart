/*
 * main.cpp
 *
 *  Created on: May 26, 2011
 *      Author: fabiorjvieira
 */

// g++ combinatoria.cpp -o combinatoria.bin
// ./combinatoria 4 AAA BBB CCC DDD EEE FFF

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fstream>

namespace Support
{
   static const unsigned int MAX_READS = 1024*256;
   static const unsigned int INVALID_IDENTIFICATION = -1;
   static const char END_OF_INFORMATION[2] = ";";
   static char SEPARATOR[2] = ",";
   static const char END_OF_LINE[2] = "\0";
   static const char COMMENTARY[2] = "#";
}

unsigned int p;
std::vector <double> pesosEntrada;
std::vector <double> pesosSaida;
std::vector <std::vector <double> > * entrada;
std::vector <std::vector <bool> > * saida;

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

void combine(std::vector <unsigned int> options, unsigned int j)
{
	unsigned int i = 0;
	options.erase(options.begin()+j);
	std::cout << options.at(0) << std::endl;
	while (options.size() > i + 1)
	{
		for (unsigned int k = 0; k < options.size(); k++) std::cout << options[k] << ",";
		std::cout << std::endl;
		combine(options,i);
		i++;
	}
}

void combinatoria(unsigned int i, std::vector <bool> vizinhos, std::vector <bool> candidatoSaida,double candidatoPeso)
{
	bool coincidencia=false;
	unsigned int k;

	for (unsigned j = i; j < entrada->size(); j++)
		vizinhos.at(j) = vizinhos.at(j) or (entrada->at(i).at(j) > 0);
	candidatoPeso+=pesosEntrada.at(i);
	candidatoSaida.at(i)=true;

	while (++i < entrada->size())
   		if (not vizinhos.at(i))
   			combinatoria(i,vizinhos,candidatoSaida,candidatoPeso);

	for (k=0; k < saida->size() and not coincidencia; k++)
	{
		coincidencia = true;
		for (unsigned int j=0; j < candidatoSaida.size(); j++)
		{
			if (not saida->at(k).at(j) and candidatoSaida.at(j))
			{
				coincidencia=false;
				break;
			}
		}
	}
//	std::cout << std::endl;
//	for (unsigned int l=0; l < candidatoSaida.size(); l++)
//		if (candidatoSaida.at(l))
//			std::cout << l << " ";
//	std::cout << "- " << candidatoPeso << std::endl;
	if (not coincidencia)
	{
		saida->push_back(candidatoSaida);
		pesosSaida.push_back(candidatoPeso);
	}
}

void permutation(std::vector < char > set, std::vector < char > combination)
{
	std::vector < char > aux;
	if (set.size() == 0)
	{
		for (unsigned int j = 0; j < combination.size(); j++)
			std::cout << combination[j];
		std::cout << std::endl;
	}
	else
		for (unsigned int j = 0; j < set.size(); j++)
		{
			combination.push_back(set[j]);
			aux = set;
			aux.erase(aux.begin()+j);
			permutation(aux, combination);
			combination.pop_back();
		}
}

void doit(unsigned int * p, unsigned int N)
{
	for (unsigned int n = 1; n <= N; n++)
		std::cout << p[n] << " ";
	std::cout << std::endl;
}

void exch(unsigned int i, unsigned int j, unsigned int * p)
{
	int t = p[i];
	p[i] = p[j];
	p[j] = t;
}

void permutation(unsigned int N)
{
	unsigned int p[N+1], c[N+1];
	for (unsigned int n = 1; n <= N; n++)
	{
		p[n] = n;
		c[n] = 1;
	}
	doit(p,N);
	for (unsigned int n = 1; n <= N;)
	{
		if (c[n] < n)
		{
			exch((n % 2) ? 1 : c[n], n, p);
			c[n]++;
			n=1;
			doit(p,N);
		}
		else c[n++] = 1;
	}
}

void permutation(unsigned int N, unsigned int * p)
{
	unsigned int c;
	if (N == 1) doit(p,4);
	else
		for (c = 1; c <= N; c++)
		{
			permutation(N-1, p);
			exch((N % 2) == 1 ? 1 : c, N, p);
		}
}

void bla(unsigned int N)
{
	unsigned int p[N+1];
	for (unsigned int n = 1; n <= N; n++) p[n] = n;
	permutation(N, p);
}



int main(int nargs, char ** args)
{
	bla(4);
	permutation(4);

	std::vector < char > set, combination;
	set.push_back('A');set.push_back('B');set.push_back('C');set.push_back('D');set.push_back('E');
	permutation(set, combination);

	/*
	std::vector <unsigned int> vt(5);
	for (unsigned int i = 0; i < vt.size(); i++)
		vt[i] = i;
	combine(vt,0);
	*/

//	saida = new std::vector <std::vector <bool> > ();
//	entrada = readCSVfile(args[1]);
//	pesosEntrada = readCSVfile(args[2])->at(0);
//	std::vector <bool> vizinhos(entrada->size(), false);
//	std::vector <bool> candidatoSaida(entrada->size(), false);
//	saida->reserve(1000000);
//
//	std::cout << "# total:" << entrada->size();
//	for (unsigned int i = 0; i < entrada->size(); i++)
//	{
//		std::cout << " " << i;
//		std::cout.flush();
//		combinatoria(i,vizinhos,candidatoSaida,0);
//	}
//
//	std::cout << std::endl;
//	for (unsigned int k=0; k < saida->size() ; k++)
//	{
//		for (unsigned int l=0; l < saida->at(k).size(); l++)
//			if (saida->at(k).at(l))
//				std::cout << l << " ";
//		std::cout << "- " << pesosSaida.at(k) << std::endl;
//	}
//
//	delete entrada;
//	delete saida;
	return 0;
}
