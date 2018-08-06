#ifndef HEADER_4D7F42FDD0838B93
#define HEADER_4D7F42FDD0838B93

#pragma once
#ifndef _GENOME_CPN_H
#define _GENOME_CPN_H

#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <fstream>


#include "ChrCopyNumber.h"
#include "SVfinder.h"
#include "myFunc.h"

class GenomeCopyNumber
{
public:
	GenomeCopyNumber(void);
	~GenomeCopyNumber(void);

	int readCGprofile(std::string const& inFile);
	void readGemMappabilityFile(std::string const& inFile);

	void printCGprofile(std::string const& outFile);
	void printCGprofile(std::string const& chr, std::ofstream & file);

	int findIndex (std::string const& chr);
	void fillCGprofile(std::string const& chrFolder);

	int getWindowSize(void);

	void setStep(int step);
    void readChrInfo( std::string const& chrLenFileName, int window, int step);

private:
	std::vector<ChrCopyNumber> chrCopyNumber_;
	std::map<std::string, int> chromosomesInd_;
	int windowSize_;
	int step_;
	double refGenomeSize_;
	bool isChrPrefix_;
};
#endif

#endif // header guard
