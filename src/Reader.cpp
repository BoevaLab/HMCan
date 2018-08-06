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
#define MAX_BUFFER 2048
#define BUFFER_SIZE 1000000

#include "Reader.h"

using namespace std;



Reader::Reader(Format type, int quality, bool duplicates) {
	// TODO Auto-generated constructor stub
	this->type = type;
	map_quality = quality;
	this->duplicates = duplicates;
}

Reader::~Reader() {
	// TODO Auto-generated destructor stub
}


reads Reader::Read(const char * filename, int & readCount){
	reads tags;
	switch(type){
		case BAM:
			read_bam(filename,tags);
			break;
		case SAM:
			read_sam(filename,tags);
			break;
		case BED:
			read_bed(filename,tags);
			break;

	}

	for(reads::iterator it=tags.begin(); it!=tags.end();++it) {
			sort((*it).second.begin(),(*it).second.end(),compare);
    }
    readCount = getNumberOfTags(tags);
	return tags;
}

int Reader::getNumberOfTags(reads&tags) {
    int counts = 0;
    reads::iterator chr_it;
	for(chr_it = tags.begin();chr_it!=tags.end();++chr_it){
        counts+=(*chr_it).second.size();
	}
	return counts;
}

reads Reader::ReadFromMultipeFiles (std::vector<std::string>& files){
	reads tags;

	switch(type){
			case BAM:
				for (unsigned int i=0;i<files.size();i++)
					read_bam(files[i].c_str(),tags);
				break;
			case SAM:
				for (unsigned int i=0;i<files.size();i++)
					read_sam(files[i].c_str(),tags);
				break;
			case BED:
				for (unsigned int i=0;i<files.size();i++)
					read_bed(files[i].c_str(),tags);
				break;

		}

		for(reads::iterator it=tags.begin(); it!=tags.end();++it)
				sort((*it).second.begin(),(*it).second.end(),compare);



	return tags;
}
char* Reader::getLine(char* buffer, int buffer_size, FILE* stream, std::string& line)
{
  while (fgets(buffer, buffer_size, stream) != NULL) {
	size_t len = strlen(buffer);
	if (buffer[len-1] == '\n') {
	  return buffer;
	}
	line = buffer;
	for (;;) {
	  if (line.at(line.length() - 1) == '\n') {
		break;
	  }
	  if (fgets(buffer, MAX_BUFFER, stream) != NULL) {
		line.append(buffer);
	  } else {
		break;
	  }
	}
	return (char*)line.c_str();
  }
  return NULL;
}
void Reader::read_bam(const char * filename,reads& tags){
	string command;
	char buffer[MAX_BUFFER];
	FILE *stream;
	string line;
	char* line_buffer;
	command = "samtools view "+string(filename);
	stream = popen(command.c_str(), "r"); //samtools should be installed in the system
	while ((line_buffer = getLine(buffer, MAX_BUFFER, stream, line)) != NULL) {
		  process_sam_line(line_buffer);
    }
    pclose(stream);
	 if (duplicates)
		 remove_duplicates();
	 merge_strands(tags);
}

void Reader::process_sam_line(const char* line_buffer) {
    if (!*line_buffer) {
        return;
    }
    string chr1;
    int orient1, start, fragmentLength, mquality,readLength;
    if (getSAMinfo(line_buffer,chr1,orient1,start, fragmentLength,mquality,readLength)) {
		DNA_fragment read;
		if (mquality > map_quality){
			read.strand = orient1;
			if (read.strand == 0){ //forward read
                read.start = start-1; //adjust to zero based counting
                if (fragmentLength>0) {
                    read.fragmentLength=fragmentLength;
                }else {
                    read.fragmentLength = 0; //as the mate is not mapped where expected (i.e. on the right)
                }
				tags_forward[chr1].push_back(read);
			}
			else if(read.strand == 1){ //reverse read
                read.start = start-1+readLength; //adjust to zero based counting
				if (fragmentLength<0) {
                    read.fragmentLength=-fragmentLength;
                }else {
                    read.fragmentLength = 0; //as the mate is not mapped where expected (i.e. on the left)
                }
				tags_reverse[chr1].push_back(read);
			}
		}
    }
}
unsigned int Reader::split(char* str_ori, char delim, char* elems[])
{
  const char* str = str_ori;
  unsigned int last_jj = 0;
  unsigned int jj = 0;
  unsigned int elem_cnt = 0;
  char c;
  for (; c = *str++; ++jj) {
	if (c == delim) {
	  str_ori[jj] = 0;
	  elems[elem_cnt++] = &str_ori[last_jj];
	  last_jj = jj+1;
	}
  }
  elems[elem_cnt++] = &str_ori[last_jj];
  return elem_cnt;
}

bool Reader::getSAMinfo(const char* line, std::string &chr1, int& orient1, int &start, int &insert_size, int & mquality, int & readLength) {

    if (!*line) {
	  return false;
	}
	if (line[0] == '@') {
	  return false;
	}

	char* strs[32];
	unsigned int strs_cnt = split((char*)line, '\t', strs);
	if (strs_cnt < 7)
        return false;

    chr1 = strs[2];
    if (chr1.compare("*")==0) {
        return false;
	}

	unsigned int mask = atoi(strs[1]);
    if (mask & 0x0010) {
	  orient1 = 1;
	} else {
	  orient1 = 0;
	}

    start = atoi(strs[3]);
    insert_size = atoi(strs[8]);
    mquality=atoi(strs[4]);
    readLength=strlen(strs[9]);

    return true;
}

void Reader::read_sam(const char * filename,reads& tags){

	ifstream infile;
	string line;
	vector<string> buf(BUFFER_SIZE);
	int i;

	infile.open(filename);
	if (!infile){
		cerr<<"Error: "<<filename<<" can not be found"<<endl;
		exit(1);
	}
	while(!infile.eof()){


		for (i=0;i<BUFFER_SIZE;i++){
			getline(infile,line);
			buf[i] = line;
			if(infile.eof())
				break;
		}
		process_sam_line(buf,i);

	}

	if (duplicates)
		remove_duplicates();
	merge_strands(tags);

}

inline void Reader::process_sam_line(vector<string>& line, int count){

	for (int i=0;i<count;i++){
		vector<string> tokens;
		DNA_fragment read;
		string chr;
		istringstream stream(line[i]);
		copy(istream_iterator<string>(stream),
			 istream_iterator<string>(),
			 back_inserter<vector<string> >(tokens)); // split the line

		if (tokens[0][0] == '@')
			continue;
		if (atoi(tokens[4].c_str()) > map_quality){
			bitset <16> flags(atoi(tokens[1].c_str()));
			read.strand = int(flags[4]);
			chr = tokens[2];
			int fragmentLength=atoi(tokens[8].c_str());
			if (read.strand == 0){ //forward read
                read.start = atoi(tokens[3].c_str())-1; //adjust to zero based counting
                if (fragmentLength>0) {
                    read.fragmentLength=fragmentLength;
                }else {
                    read.fragmentLength = 0; //as the mate is not mapped where expected (i.e. on the right)
                }
				tags_forward[chr].push_back(read);
			}
			else if(read.strand == 1){ //reverse read
                int readLength=tokens[9].length();
                read.start = atoi(tokens[3].c_str())-1+readLength; //adjust to zero based counting
				if (fragmentLength<0) {
                    read.fragmentLength=-fragmentLength;
                }else {
                    read.fragmentLength = 0; //as the mate is not mapped where expected (i.e. on the left)
                }
				tags_reverse[chr].push_back(read);
			}
		}
	}

}

void Reader::read_bed(const char * filename,reads& tags){
	ifstream infile;
		string line;
		map<string,set<long> >::iterator chr_it;
		int i;
		vector<string> buf(BUFFER_SIZE);
		infile.open(filename);
		if (!infile){
			cerr<<"Error: "<<filename<<" can not be found"<<endl;
			exit(1);
		}




	while(!infile.eof()){


			for (i=0;i<BUFFER_SIZE;i++){
				getline(infile,line);
				buf[i] = line;
				if(infile.eof())
					break;
			}
			process_bed_line(buf,i);

		}
		if (duplicates)
			remove_duplicates();
		merge_strands(tags);
}
inline void Reader::process_bed_line(vector<string>& line, int count){

	for (int i=0;i<count;i++){
		vector<string> tokens;
		DNA_fragment read;
		string chr;
		istringstream stream(line[i]);
		copy(istream_iterator<string>(stream),
			 istream_iterator<string>(),
			 back_inserter<vector<string> >(tokens)); // split the line
		if (tokens[0].compare("track")==0)
			return; //skip the track header
		read.strand = tokens[5].compare("+")==0 ? 0:1;
		chr = tokens[0];
		read.start = read.strand == 0 ? atoi(tokens[1].c_str()):atoi(tokens[2].c_str())-1;
		if (read.strand == 0){
				tags_forward[chr].push_back(read);
		}

		else if (read.strand == 1){
			tags_reverse[chr].push_back(read);
		}

	}


}

void Reader::remove_duplicates(){
	reads::iterator chr_it;
	for(chr_it = tags_forward.begin();chr_it!=tags_forward.end();++chr_it){
		sort((*chr_it).second.begin(),(*chr_it).second.end(),compare);
		(*chr_it).second.erase(unique((*chr_it).second.begin(),(*chr_it).second.end(),equivelant),
				(*chr_it).second.end());
	}

	for(chr_it = tags_reverse.begin();chr_it!=tags_reverse.end();++chr_it){
			sort((*chr_it).second.begin(),(*chr_it).second.end(),compare);
			(*chr_it).second.erase(unique((*chr_it).second.begin(),(*chr_it).second.end(),equivelant),
					(*chr_it).second.end());
	}
}


void Reader:: merge_strands(reads& tags){

	reads::iterator chr_it;

	for(chr_it=tags_forward.begin();chr_it!=tags_forward.end();++chr_it){
		for(unsigned int i=0;i<(*chr_it).second.size();i++)
			tags[(*chr_it).first].push_back((*chr_it).second[i]);
		(*chr_it).second.clear();
	}

	for(chr_it=tags_reverse.begin();chr_it!=tags_reverse.end();++chr_it){
			for(unsigned int i=0;i<(*chr_it).second.size();i++)
				tags[(*chr_it).first].push_back((*chr_it).second[i]);
			(*chr_it).second.clear();
	}
}


vector<BedEntry> Reader::Read_blacklist(const char * filename){
	BedEntry temp;
	vector<BedEntry> blacklist;
	ifstream file;

	if (filename=="")
        return blacklist;
	file.open(filename);
	if(!file){
		cerr<<"Can not open blacklist file. Exiting..."<<endl;
	}

	do{
		file>>temp.chr;
		file>>temp.start;
		file>>temp.end;
		if (file.eof())
			break;
		blacklist.push_back(temp);
	} while(1);

	return blacklist;
}
