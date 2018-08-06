#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#define NA -1
#include "types.h"


bool equivelant(DNA_fragment tag1, DNA_fragment tag2);
bool compare(DNA_fragment tag1, DNA_fragment tag2);
bool compare_2d_vector(std::vector<float> , std::vector<float>);
std::vector <std::string > ReadFasta(const char * file_name);
std::string Read_chr(const char * file_name);
float get_max_element(float * values, int, int);
//float get_median(const std::vector<float> &);
float calculate_area(float*, int ,int);
std::string long2str(long);
std::string float2str(float);
bool BedCompare(BedEntry a, BedEntry b);
bool BedStartCompare(BedEntry a, BedEntry b);
int iround (float n);
bool is_file_exist(const char *fileName);

#endif /* UTILS_H_ */
