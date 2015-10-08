/*
 * main.hpp
 *
 *  Created on: Jun 5, 2010
 *      Author: rocha
 */

#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <stdlib.h>
#include <algorithm>

#include <iostream>
#include <sstream>
#include <fstream>

#include <cmath>
#include <vector>

#include "./Statistics/studenttdistr.h"

namespace Support
{
   static const unsigned int MAX_READS = 1024*256;
   static const unsigned int INVALID_IDENTIFICATION = -1;
   static const char END_OF_INFORMATION[2] = ";";
   static const char SEPARATOR[2] = ",";
   static const char END_OF_LINE[2] = "\0";
   static const char COMMENTARY[2] = "#";
   static const double CONFIDENCE_FACTOR = 0.95;
}

#endif /* MAIN_HPP_ */
