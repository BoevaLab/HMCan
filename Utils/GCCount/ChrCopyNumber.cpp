/*************************************************************************
Copyright (c) 2010-2011, Valentina BOEVA.

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


#include "ChrCopyNumber.h"

using namespace std ;

ChrCopyNumber::ChrCopyNumber(std::string const& chrName) {
	chromosome_ = chrName;
}

ChrCopyNumber::ChrCopyNumber(void) {
}

ChrCopyNumber::ChrCopyNumber(int windowSize, int chrLength, std::string const& chrName) {
	windowSize_ = windowSize;
	chrLength_ = chrLength;
	chromosome_ = chrName;
	step_=windowSize;
	if (windowSize ==0) {
        cerr << "Error: windowSize is set to Zero\n";
        exit(-1);
	}
	length_ = chrLength/windowSize+1;
	coordinates_ = vector<int>(length_);
	for (int i = 0; i<length_; i++) {
		coordinates_[i] = i*windowSize;
	}

}

ChrCopyNumber::ChrCopyNumber(int windowSize, int chrLength, std::string const& chrName, int step) {
	windowSize_ = windowSize;
	step_=step;
	chrLength_ = chrLength;
	chromosome_ = chrName;

	if (windowSize ==0) {
        cerr << "Error: windowSize is set to Zero\n";
        exit(-1);
	}
	length_ = chrLength/step+1;
	coordinates_ = vector<int>(length_);
	if (step<windowSize) {
        ends_ = vector<int>(length_,0);
        for (int i = 0; i<length_; i++) {
            ends_[i] = i*step+windowSize_-1;
        }
	}
	for (int i = 0; i<length_; i++) {
		coordinates_[i] = i*step;
	}
}


int ChrCopyNumber::getCoordinateAtBin(int i) {
	return coordinates_[i];
}
int ChrCopyNumber::getEndAtBin(int i) {
	if ((int)ends_.size()>i)
        return ends_[i];
    else
        return coordinates_[i]+windowSize_-1;
}

void ChrCopyNumber::setNotNprofileAt(int i, float value) {
    notNprofile_[i] = value;
}

void ChrCopyNumber::setMappabilityProfileAt(int i, float value) {
    mappabilityProfile_[i] = value;
}



void ChrCopyNumber::addToCoordinates(int i) {
	coordinates_.push_back(i);
}

void ChrCopyNumber::addToEnds(int i) {
	ends_.push_back(i);
}


void ChrCopyNumber::setWindowSize(int windowSize) {
	windowSize_ = windowSize;
}


void ChrCopyNumber::setStep(int step) {
	step_ = step;
}

void ChrCopyNumber::setVectorLength(int length){
	length_ = length;
}
void ChrCopyNumber::setChrLength(int chrLength) {
	chrLength_ = chrLength;
}


float ChrCopyNumber::getCGprofileAt(int i) {
	return GCprofile_[i];
}
float ChrCopyNumber::getNotNprofileAt(int i) {
	return notNprofile_[i];
}

float ChrCopyNumber::getMappabilityProfileAt(int i) {
	return mappabilityProfile_[i];
}

int ChrCopyNumber::getLength() {
	return length_;
}
int ChrCopyNumber::getChrLength() {
	return chrLength_;
}

int ChrCopyNumber::getMappabilityLength() {
	return mappabilityProfile_.size();
}


std::string ChrCopyNumber::getChromosome() {
	return chromosome_;
}

void ChrCopyNumber::clearCGcontent () {
	GCprofile_.clear ();
}
void ChrCopyNumber::clearNonNpercent () {
	notNprofile_.clear ();
}

void ChrCopyNumber::clearMappabilityProfile () {
	mappabilityProfile_.clear ();
}

void ChrCopyNumber::addToCGcontent (float valueToAdd) {
	GCprofile_.push_back(valueToAdd);
}
void ChrCopyNumber::addToNonNpercent (float valueToAdd) {
	notNprofile_.push_back(valueToAdd);
}

void ChrCopyNumber::addToMappabilityProfile (float valueToAdd) {
	mappabilityProfile_.push_back(valueToAdd);
}

void ChrCopyNumber::createMappabilityProfile() {
    for (int i=0; i<length_; i++)
        mappabilityProfile_.push_back(0);
}

void ChrCopyNumber::checkOrCreateNotNprofileWithZeros() {
    int sizeOfNonN = notNprofile_.size();
    if (sizeOfNonN<length_) {
        notNprofile_ = vector <float> (length_);
    }
    for (int i = 0; i< length_; i++)
        notNprofile_[i]=0;
}

void ChrCopyNumber::fillCGprofile(std::string const& chrFolder) {
	GCprofile_ = vector <float> (length_);
	notNprofile_ = vector <float> (length_);
	ifstream file;
	string filename = chromosome_;
	string possibleFilenames[] = {filename,filename+".fa",filename+".fasta","chr"+filename+".fa","chr"+filename+".fasta"};
	for (int i = 0; i < 5; i++) {

		string myFilename = possibleFilenames[i];
		string myFullPath = pathAppend(chrFolder,myFilename);
		file.open(myFullPath.c_str());

		if(!file.is_open())
			file.clear();
		else
			i = 6;
	}


    if (!file.is_open())	 {
		//	throw ("Unable to open fasta file for chr "+chromosome_+" in "+chrFolder+"\n");
        cerr << "Unable to open fasta file for chr "+chromosome_+" in folder "+chrFolder+"\n\nPlease note, "<< chrFolder << " should be a folder, not a file!\n\n";
        exit (-1);
	}
	//string myString;
	char letter; // here we read the chromosome letter by letter. Can do it better!
	file >> letter;
	int count = 0;
	int countCG = 0;
	int countN = 0;
	string line;
    string text = "";

	if (letter == '>')
		getline (file,line);
	else {
		count = 1;
		countCG = isCG(letter);
		countN = isN(letter);
		text.push_back(letter);
	}

	if (ends_.size()==0) { //all windows have equal length => can use the same windowsize for all windows
		for (int i = 0; i<length_; i++) {
			if (file.eof()) {
				GCprofile_[i] = NA;
				//cout << "End-of-file reached.." << endl;
			}
			while((!file.eof()) && (count < windowSize_)) {
				file>>letter;
				countCG += isCG(letter);
				countN += isN(letter);
				count ++;
			}
			notNprofile_[i] = float(count-countN)/count;
			if (count == countN)
				GCprofile_[i] = NA;
			else
				GCprofile_[i] = float(countCG)/(count-countN);
			//reset
			countCG = 0;
			countN = 0;
			count = 0;
		}
	} else {
		int start, end;
		for (int i = 0; i<length_; i++) {
			if (file.eof()) {
				GCprofile_[i] = NA;
				//cout << "End-of-file reached.." << endl;
			}
			start = coordinates_[i];
			end = ends_[i];
			while((!file.eof()) && (count < start)) {
				file>>letter;
				count ++;
			}
			while((!file.eof()) && (count <= end)) {
				file>>letter;
				text.push_back(letter);
				count ++;
			}
			notNprofile_[i] = 1;
			//notNprofile_[i] = float(end-start+1-countN)/(end-start+1);
			/*if (end-start+1 == countN)
				GCprofile_[i] = NA;
			else */

			countCG =0;
			countN = 0;
            for(int j = 0; j < (int)text.length(); j++) //++j????
                if (text[j] == 'C' || text[j] == 'G' || text[j] == 'c' || text[j] == 'g')
                    countCG++;
                else if (text[j] == 'N')
                    countN++;
            if (end-start+1-countN>0)
                GCprofile_[i] = float(countCG)/(end-start+1-countN);
            else
                GCprofile_[i] = NA;
			//reset
			if (i+1<length_) {
			    int nextStart = coordinates_[i+1];
			    if (nextStart<=end) {
                    //count = end;
                    //and delete prefix in text;
                    int howMuchToDelete = nextStart - start;
                    text = text.substr(howMuchToDelete);
			    } else {
                    text = "";
			    }
			}

		}
	}



	file.close();

}


int ChrCopyNumber::getEndsSize() {
	return ends_.size();
}

ChrCopyNumber::~ChrCopyNumber(void)
{
	coordinates_.clear();
	length_ = 0;
}
