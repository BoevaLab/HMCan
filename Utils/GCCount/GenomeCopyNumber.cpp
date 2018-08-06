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


#include "GenomeCopyNumber.h"

using namespace std ;

GenomeCopyNumber::GenomeCopyNumber(void)
{
	step_=NA;
}

void GenomeCopyNumber::readChrInfo( std::string const& chrLenFileName, int windowSize, int step) {
    if (step == NA)
        step = windowSize;
    if (step <= 0 || step > windowSize) {
        cerr << "step  should be a positive interger value less than or equal to the window size\n";
        exit(-1);
    }
    step_=step;
    //reading the file with genome information
	std::vector<std::string> names;
	std::vector<int> lengths;
	isChrPrefix_=readFileWithGenomeInfo(chrLenFileName, names, lengths);
	refGenomeSize_ = sum(lengths);
	for (int i = 0; i < (int) names.size(); i++) {
		ChrCopyNumber chrCopyNumber(windowSize, lengths[i],names[i], step);
        //cout << names[i] << "\t" << i << "\n";
		chromosomesInd_.insert(pair<string, int> (names[i],i));
		chrCopyNumber_.push_back(chrCopyNumber);
	}
}

int GenomeCopyNumber::getWindowSize(void) {
	return windowSize_;
}

void GenomeCopyNumber::setStep(int step) {
    step_=step;
}


int GenomeCopyNumber::findIndex (std::string const& chr) {
	if (chromosomesInd_.find(chr) == chromosomesInd_.end()) {return NA;}
	return chromosomesInd_.find(chr)->second;
}

void GenomeCopyNumber::readGemMappabilityFile(std::string const& inFile) {
	ifstream file (inFile.c_str());
	string line;
	string chromInLine;
    string::size_type pos = 0;
    string currentChr = "";
    int index = NA;
    int count = 0;
    int uniqueCount = 0;
    int localWindowSize = windowSize_;
    float ratio;
    int positionInd = 0;
    int startPos = 0;
    int endPos = 0;
    int lastEnd = -1;
    bool chrstart = 0;
    string text = "";
	if (file.is_open())	{
	    cout << "..Reading "<< inFile << "\n";
		while (! file.eof() )	{
			getline (file,line);
			if (! line.length()) continue;
			pos = 0;
			if ( ( pos = line.find("~chr", pos)) != string::npos ){
			   chromInLine = line.substr(4);
			   chrstart = 1;
			}
			else if ( ( line.find("~~") == string::npos) && ( line.substr(0,1)=="~") ){
			   chromInLine = line.substr(1);
			   chrstart = 1;
			   if (( pos = line.find(" ")) != string::npos) {
                    chromInLine = line.substr(1,pos-1);
			   }
			}
            if ( chrstart ){
                //save the last info:
                if (count > endPos && index != NA) { // endPos ?
                    int howMuchToDelete = startPos-lastEnd-1;
                    text = text.substr(howMuchToDelete);
                    uniqueCount =0;
                    for(int i = 0; i < (int)text.length(); ++i)
                        if (text[i] == '!')
                            uniqueCount++;

                    ratio = float(uniqueCount)/localWindowSize;
                    chrCopyNumber_[index].setMappabilityProfileAt(positionInd, ratio);
                }
                //restore all variables for a new chromosome
                currentChr = chromInLine;
                chrstart = 0;
                cout << "..Reading mappability for chromosome " << currentChr << "\n";
                index = findIndex(currentChr);

                positionInd = 0;
				count = 0;
				uniqueCount = 0;
                text = "";
                lastEnd = -1;

                //cout <<  "..Index for chromosome " << currentChr << ": "<< index << "\n";

				if (index == NA) {
				    cout <<  "skipping chromosome " << currentChr << "\n";
				    //return;
				    // do not return! they can be other "good" chromosomes afterwords!
				    //do nothing!!! wait for the next chromosome!

                } else {
                    startPos = chrCopyNumber_[index].getCoordinateAtBin(positionInd);
                    endPos = chrCopyNumber_[index].getEndAtBin(positionInd);
                    localWindowSize = endPos-startPos+1;
                    int maxInd = chrCopyNumber_[index].getLength()-1;
                    cout << "..Control: Last window: from " << chrCopyNumber_[index].getCoordinateAtBin(maxInd) << " to " << chrCopyNumber_[index].getEndAtBin(maxInd) <<"\n";
                    chrCopyNumber_[index].createMappabilityProfile();
                }

			} else if (index != NA) {
                count += line.length();
                text.append(line);

                if (count > endPos) {
                    int howMuchToDelete = startPos-lastEnd-1;
                    if (howMuchToDelete<0)
                        howMuchToDelete=0; //this should never happen
                    text = text.substr(howMuchToDelete);
                    lastEnd = endPos;
                    string substr_ = text.substr(0,localWindowSize);
                    uniqueCount = 0;
                    for(int i = 0; i < (int)substr_.length(); ++i)
                        if (substr_[i] == '!')
                            uniqueCount++;


                    if (positionInd+1<chrCopyNumber_[index].getLength()) {
                        int nextStart = chrCopyNumber_[index].getCoordinateAtBin(positionInd+1);
                        if (nextStart<=endPos) {
                            //count = end;
                            //and delete prefix in text;
                            int howMuchToDelete = nextStart - startPos;
                            text = text.substr(howMuchToDelete);
                        } else {
                            text = text.substr(localWindowSize); //actually this should be the same as "text.substr(howMuchToDelete);"
                        }

                    }

                    //count -= localWindowSize;

                    ratio = float(uniqueCount)/localWindowSize;
                    chrCopyNumber_[index].setMappabilityProfileAt(positionInd, ratio);
                    positionInd ++;

                    if (positionInd<chrCopyNumber_[index].getLength()) {
                        startPos = chrCopyNumber_[index].getCoordinateAtBin(positionInd);
                        endPos = chrCopyNumber_[index].getEndAtBin(positionInd);
                        localWindowSize = endPos-startPos+1;
                    } else {
                        //should not read this chromosome any more
                        index = NA;
                    }
                }
			}
		}
		file.close();
		//save the very last info (int text variable)
        if (count >= startPos && index != NA) { //endPos
            uniqueCount =0;
               for(int i = 0; i < (int)text.length(); ++i)
                   if (text[i] == '!')
                        uniqueCount++;
            ratio = float(uniqueCount)/localWindowSize; //  /max((int)text.length(),(int)localWindowSize);
            chrCopyNumber_[index].setMappabilityProfileAt(positionInd, ratio);
        }

		cout << "file " << inFile << " is read\n";
//		if (chrCopyNumber_[1].getLength() == chrCopyNumber_[1].getMappabilityLength())
//            cout << "Mappability profile has been set correctly";
	} else{
        cerr << "Error: Unable to open file "+inFile+"\n";
        exit(-1);
	}
}

int GenomeCopyNumber::readCGprofile(std::string const& inFile) {
	ifstream file (inFile.c_str());
	string line;
	int count = 0;
	int observedStep = 0;
	if (file.is_open())	{
		while (! file.eof() )	{
			getline (file,line);
			if (! line.length()) continue;
			std::vector<std::string> strs = split(line, '\t');
			if (strs.size()>=4) {
				string currentChr = strs[0];
				if (observedStep==0)
                    observedStep = ceil(strtod(strs[1].c_str(), NULL));
				float CGperc =(float)strtod(strs[2].c_str(), NULL);
				float nonNperc =(float)strtod(strs[3].c_str(), NULL);
				string::size_type pos = 0;
				if ( ( pos = currentChr.find("chr", pos)) != string::npos )
					currentChr.replace( pos, 3, "" );
				int index = findIndex(currentChr);
				if (index != NA) {
				    chrCopyNumber_[index].addToCGcontent(CGperc);
                    chrCopyNumber_[index].addToNonNpercent(nonNperc);
                    count++;
                    if (strs.size()==5) { //means that there are also mappability values in the 5th colomn
                        float MappPerc =(float)strtod(strs[4].c_str(), NULL);
                        chrCopyNumber_[index].addToMappabilityProfile(MappPerc);
                    }
                    if (observedStep!=0 && observedStep!=step_) {
                        file.close();
                        chrCopyNumber_[index].clearCGcontent();
                        chrCopyNumber_[index].clearNonNpercent();
                        chrCopyNumber_[index].clearMappabilityProfile();
                        return observedStep;

                    }
				}
			}

			strs.clear();
		}
		file.close();
		cout << "file " << inFile << " is read\n";
		if (count==0){
            cerr << "Your GC-content file "<<inFile<< " is empty or is in a wrong format\n\nPlease use chomosome sequences (option \"chrFiles\") to recreate it!\n\n";
            exit(-1);
		}
	} else {
	    cerr << "Unable to open file "+inFile+"\n";
	    exit (-1);
    }

    return observedStep;
}


void GenomeCopyNumber::printCGprofile(std::string const& outFile) {
	const char * name = outFile.c_str();
	std::ofstream file;
	file.open(name);
	map<string,int>::iterator it;
	for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
		printCGprofile((*it).first,file);
	}
	file.close();
	cout << "CG-content printed into "<<outFile <<"\n";
}


void GenomeCopyNumber::printCGprofile(std::string const& chr, std::ofstream & file) {
	string::size_type pos = 0;
	string chrNumber = chr;
	string chrPrefix = "";
	if ( ( pos = chrNumber.find("chr", pos)) != string::npos ) {
		chrNumber.replace( pos, 3, "" );
		chrPrefix = "chr";
	}
	if (isChrPrefix_){
		chrPrefix = "chr";
	}
	map<string,ChrCopyNumber>::iterator it;
	int index = findIndex(chrNumber);
	if (index == NA) {return;}
	int length = chrCopyNumber_[index].getLength();
	if (chrCopyNumber_[index].getMappabilityLength()>0)
        for (int i = 0; i< length; i++)
            file << chrPrefix<<chrNumber <<"\t"<< chrCopyNumber_[index].getCoordinateAtBin(i) <<"\t"<< chrCopyNumber_[index].getCGprofileAt(i) <<"\t"<<chrCopyNumber_[index].getNotNprofileAt(i) <<"\t"<< chrCopyNumber_[index].getMappabilityProfileAt(i) << "\n";
    else
        for (int i = 0; i< length; i++)
            file << chrPrefix<<chrNumber <<"\t"<< chrCopyNumber_[index].getCoordinateAtBin(i) <<"\t"<< chrCopyNumber_[index].getCGprofileAt(i) <<"\t"<<chrCopyNumber_[index].getNotNprofileAt(i) << "\n";
}


void GenomeCopyNumber::fillCGprofile(std::string const& chrFolder) {
	//reading the file with genome information
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		it->fillCGprofile(chrFolder);
	}
}

GenomeCopyNumber::~GenomeCopyNumber(void)
{
	chrCopyNumber_.clear();
	chromosomesInd_.clear();
}
