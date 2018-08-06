/*
 * types.h
 *
 *  Created on: Sep 7, 2014
 *      Author: haitham
 */

#ifndef TYPES_H_
#define TYPES_H_

#include<vector>
#include <map>
#include <string>
enum Format {
	BAM,
	BED,
	SAM,
  WIG
};



struct DNA_fragment{
	int start;
	bool strand; // 0 for forward and 1 for reveresd
	int fragmentLength;
};

struct BedEntry{
	std::string chr;
	int start;
	int end;
	float score;
	int max;
	int max_value;
	double p_value;

};

struct GC_vectors{
	std::vector<int> starts;
	std::vector<int> ends;
	std::vector<float> GC_content;
};


struct DNA_bin{
	int start;
	std::vector<float> density;
};

struct observation_seq{
	std::string chr;
	int start;
	std::vector<int> values;
	std::vector<double> posterior_prob;

};



typedef std::map<std::string,std::vector<float> > density;
typedef std::map<std::string,std::vector<DNA_fragment> > reads;
typedef std::map<std::string,std::vector<int> > counts ;
typedef std::vector<std::vector<float> > table_2d;
typedef std::vector<std::vector<int> > table_2d_i;
typedef std::map<std::string,std::vector<BedEntry> > bed_map;
typedef std::vector<density> density_vector;
typedef std::vector<density_vector> density_matrix;

/*those two can be optimized*/
typedef std::vector<std::vector<counts> > counts_matrix;
typedef std::vector<counts> counts_vector;
/*************************************************/


typedef std::vector<reads> reads_vector;
typedef std::vector<std::vector<reads> > reads_matrix;
#endif /* TYPES_H_ */
