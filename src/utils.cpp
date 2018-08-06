/*************************************************************************
Copyright (c) 2013, Haitham ASHOOR

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#include "utils.h"

using namespace std;

bool compare(DNA_fragment tag1, DNA_fragment tag2){
	return tag1.start<tag2.start;
}

bool equivelant(DNA_fragment tag1, DNA_fragment tag2){
	return tag1.start == tag2.start;
}
vector <string > ReadFasta(const char * file_name){

	ifstream infile;

	string str,temp;
	vector <std::string > seqs;
	infile.open(file_name);
	if (!infile)
		std::cout<<"unable to open the file"<<std::endl;


	do {
		getline(infile,temp);
		if (infile.eof())
				break;
		if (temp[0] == '>')
		{
			if (!str.empty()){

				transform(str.begin(), str.end(),str.begin(), ::toupper);
				seqs.push_back(str);
				str.clear();// = '';
			}
		}
		else if (temp[0] != '>' && !temp.empty())
		{

			str = str+temp;
			temp.clear();

		}





	}while(true);
	if (!str.empty())
			seqs.push_back(str);
	return seqs;
}


bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

string Read_chr(const char * file_name){
	string chr,line;
	vector <string> lines;
	ifstream infile;
	infile.open(file_name);

	if (!infile){
		cerr<<"Error: Can not open file "<<file_name<<" Exiting...."<<endl;
		exit(1);
	}
	chr = "";
	getline(infile,line);

	while(getline(infile,line)){
		lines.push_back(line);
	}

	for(unsigned int i=0;i<lines.size();i++)
		chr.append(lines[i]);
	return chr;
}


float get_max_element(float * values, int start, int end){
	float max = values[start];
	for (int i=1;i<end;i++)
		if(values[i]>max)
			max = values[i];
	return max;
}



float calculate_area(float *values, int start,int end){ // calculate area under curve using trapozoidal rule
	float norm_factor = 1/(2*(end-start));
	float area = norm_factor*(values[start]+values[end]);
	for(int i=start+1;i<end;i++)
		area+=2*norm_factor*values[i];
	return area;
}


string long2str(long number){
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

string float2str(float number){
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}


bool compare_2d_vector(vector<float> a,vector<float> b){
	return a[0]<b[0];
}


bool BedCompare(BedEntry a, BedEntry b){
	return a.score>b.score;
}

bool BedStartCompare(BedEntry a, BedEntry b){
	return a.start<b.start;
}



int iround(float n){
	n = n*100;
	if (int(n) - n >0.5)
		return (int) n +1;
	else
		return (int) n;



}
