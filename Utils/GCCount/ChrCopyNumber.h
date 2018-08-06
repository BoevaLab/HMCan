#ifndef HEADER_6AE83A503DDCA686
#define HEADER_6AE83A503DDCA686

#pragma once
#ifndef _CHR_CPN_H
#define _CHR_CPN_H

#include <vector>
#include <map>
#include <stdlib.h>
#include <iostream>

#include "myFunc.h"
//#include "SVfinder.h"

class ChrCopyNumber
{
public:
    ChrCopyNumber(void);
	ChrCopyNumber(std::string const& chrName);
	ChrCopyNumber(int windowSize, int chrLength, std::string const& chrName);
	ChrCopyNumber(int windowSize, int chrLength, std::string const& chrName, int step);
	~ChrCopyNumber(void);

	float		getValueAt(int i);
	int			getCoordinateAtBin(int i);
	int			getEndAtBin(int i);
		int			getLength();
    int         getChrLength();
	std::string getChromosome();

	float getCGprofileAt(int i);
	float getMappabilityProfileAt(int i);
	float getNotNprofileAt(int i);


	int getEndsSize();
	int	getMappabilityLength();

	void fillCGprofile(std::string const& chrFolder);

	void addToReadCount(float);
	void addToCoordinates(int);
	void addToEnds(int i);
	void addToCGcontent (float valueToAdd);
	void addToNonNpercent (float valueToAdd);
    void addToMappabilityProfile(float valueToAdd);

    void clearCGcontent () ;
    void clearNonNpercent () ;
    void clearMappabilityProfile ();

	void setWindowSize(int);
	void setVectorLength(int);
	void setChrLength(int);
	void setNotNprofileAt(int i, float value);
	void setMappabilityProfileAt(int i, float value);
	void setStep(int step);
	void setRCountToZeroForNNNN();
    void createMappabilityProfile();
    void checkOrCreateNotNprofileWithZeros();
	float getLevelAt (int unsigned i, int ploidy);


private:
    int ploidy_;
    float normalContamination_;
	int chrLength_;
	int length_;
	int windowSize_;
	int step_;
	bool isMedianCalculated_;
	bool isSmoothed_; //medianProfileHasBeenSmoothed
	std::string chromosome_;
	std::vector <int> coordinates_;
	std::vector <int> ends_;
	std::vector <int> bpfinal_;
	std::vector <int> fragmentNotNA_lengths_;
	std::vector <int> fragment_lengths_;
	std::vector <float> GCprofile_;//CG-content in a window
	std::vector <float> notNprofile_;//percentage of not 'N' in a window
	std::vector <float> mappabilityProfile_;//percentage of mappable positions in a window
};
#endif

#endif // header guard
