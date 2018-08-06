#ifndef PARSER_H_
#define PARSER_H_
#include <string>
#include "utils.h"
class Parser {
public:
	Parser();
	virtual ~Parser();
	void parse(char *);
	void print();
	int min_length, median_length,max_length,bin_size,large_bin_size,merge_dist,max_iter, GC_merge_dist;
	float iter_threshold,final_threshold,pvalue_threshold, posterior_threshold;
	bool wig,posterior,bedgraph,callpeaks,pairedEnds,removeDuplicates,isCHIPseq,calculateEmpiricalPvalue,ifKeepChIPdensityIntact,CNAnormalization;
	std::string path,GC_index,blacklist,genomeLengthFile;
	Format t;

};

#endif /* PARSER_H_ */
