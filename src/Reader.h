#ifndef READER_H_
#define READER_H_

#include <vector>
#include <string>
#include <map>
#include <bitset>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <set>
#include<cstdlib>
#include "utils.h"
#include "types.h"
#include <string.h>

class Reader {
public:
	Reader(Format, int, bool);
	virtual ~Reader();
	reads Read(const char *, int &);
	reads ReadFromMultipeFiles (std::vector<std::string>&);
	std::vector<BedEntry> Read_blacklist(const char *);
	bool getSAMinfo(const char* line, std::string &chr1, int& orient1, int &start, int &insert_size, int & mquality, int & readLength);
    char* getLine(char* buffer, int buffer_size, FILE* stream, std::string& line);
    unsigned int split(char* str_ori, char delim, char* elems[]);

private:
	Format type;
	int map_quality;
	bool duplicates;
	int getNumberOfTags(reads&) ;
	reads tags_forward, tags_reverse;
	void read_bam(const char * ,reads&);
	void read_sam(const char * ,reads&);
	void read_bed(const char * ,reads&);
	inline void process_sam_line(std::vector<std::string>&,int);
	inline void process_bed_line(std::vector<std::string>&,int);
	void remove_duplicates();
	void merge_strands(reads&);
	void process_sam_line(const char* line_buffer) ;

};

#endif /* READER_H_ */
