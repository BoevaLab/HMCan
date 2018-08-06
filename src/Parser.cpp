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

#include "Parser.h"
#include <fstream>
#include <iostream>
#include<algorithm>
using namespace std;
#include <cstdlib>

Parser::Parser() {
	GC_index="";
	path ="";
	blacklist = "";
	genomeLengthFile="";
	min_length=0;
	median_length=0;
	max_length=0;
	large_bin_size=0;
	bin_size=0;
	pvalue_threshold =-1;
	merge_dist = -1;
	iter_threshold = -1;
	final_threshold = -1;
	max_iter = 0;
	posterior_threshold =-1;
	wig = false;
	posterior = false;
	t = BAM;
	callpeaks=true;
    pairedEnds=false;
    removeDuplicates=false;
    isCHIPseq=true;
    GC_merge_dist=1000;
    calculateEmpiricalPvalue=false; //TODO: may want to change it to true
    bedgraph=false;
    ifKeepChIPdensityIntact=true;
    CNAnormalization=true;

}

Parser::~Parser() {
	// TODO Auto-generated destructor stub
}

void Parser::parse(char * config_file){

	ifstream conf;
	conf.open(config_file);
	string item,tempStr;

	if (!conf){
		cerr<<"Can not open the file "<<config_file<<endl;
		exit(1);
	}
		//geting the configuration
		while(!conf.eof()){

			conf>>item>>tempStr;
			if (item.compare("GCIndex")==0)
				GC_index = tempStr;
			else if (item.compare("genomePath")==0)
				path = tempStr;
			else if (item.compare("minLength")==0)
				min_length = atoi(tempStr.c_str());
			else if(item.compare("medLength")==0)
				median_length = atoi(tempStr.c_str());
			else if (item.compare("maxLength")==0)
				max_length= atoi(tempStr.c_str());
			else if (item.compare("smallBinLength")==0)
				bin_size = atoi(tempStr.c_str());
			else if (item.compare("largeBinLength")==0)
				large_bin_size = atoi(tempStr.c_str());
            else if (item.compare("genomeLengthFile")==0)
				genomeLengthFile = tempStr;
			else if (item.compare("pvalueThreshold")==0)
				pvalue_threshold = atof(tempStr.c_str());
			else if (item.compare("GCmergeDistance")==0)
				GC_merge_dist = atoi(tempStr.c_str());
            else if (item.compare("mergeDistance")==0)
				merge_dist = atoi(tempStr.c_str());
			else if(item.compare("iterationThreshold")==0)
				iter_threshold = atof(tempStr.c_str());
			else if (item.compare("finalThreshold")==0)
				final_threshold = atof(tempStr.c_str());
			else if (item.compare("maxIter")==0)
				max_iter = atoi(tempStr.c_str());
			else if (item.compare("format")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("bed")==0)
					t = BED;
				else if (tempStr.compare("sam")==0)
					t = SAM;
				else if (tempStr.compare("bam")==0)
					t = BAM;
                else if (tempStr.compare("wig")==0)
                    t = WIG;
				else{
					cerr<<"Please provide a file in SAM or BED format."<<endl;
					exit(EXIT_FAILURE);
				}
			}
			else if (item.compare("PrintWig")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("true")==0)
					wig = true;
				else if (tempStr.compare("false")==0)
					wig = false;
				else{
					wig = false;
					cerr<<"Warning: "<<item<<" is not a valid option for PrintWig; it has to be TRUE or FALSE"<<endl;
					cerr<<"Warning: PrintWig is set to FALSE"<<endl;
				}
			}
            else if (item.compare("calculateEmpiricalPvalue")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("true")==0) {
					calculateEmpiricalPvalue = true;
                    cout<<"Warning: calculateEmpiricalPvalue was set to true; The calculations will take twice more time!"<<endl;

				}
				else if (tempStr.compare("false")==0)
					calculateEmpiricalPvalue = false;
				else{
					calculateEmpiricalPvalue = false; //TODO: modify to true if the option works well in the tests
					cerr<<"Warning: "<<item<<" is not a valid option for calculateEmpiricalPvalue; it has to be TRUE or FALSE"<<endl;
					cerr<<"Warning: calculateEmpiricalPvalue is set to FALSE"<<endl;
				}
			}
			else if (item.compare("keepChIPdensityIntact")==0) {
                transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("true")==0)
					ifKeepChIPdensityIntact = true;
				else if (tempStr.compare("false")==0)
					ifKeepChIPdensityIntact = false;
				else{
					ifKeepChIPdensityIntact = true;
					cerr<<"Warning: "<<item<<" is not a valid option for keepChIPdensityIntact; it has to be TRUE or FALSE"<<endl;
					cerr<<"Warning: keepChIPdensityIntact is set to TRUE"<<endl;
				}
			}
            else if (item.compare("CNAnormalization")==0) {
                transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("true")==0)
					CNAnormalization = true;
				else if (tempStr.compare("false")==0)
					CNAnormalization = false;
				else{
					CNAnormalization = true;
					cerr<<"Warning: "<<item<<" is not a valid option for CNAnormalization; it has to be TRUE or FALSE"<<endl;
					cerr<<"Warning: CNAnormalization is set to TRUE"<<endl;
				}
			}
			else if (item.compare("RemoveDuplicates")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("true")==0)
					removeDuplicates = true;
				else if (tempStr.compare("false")==0)
					removeDuplicates = false;
				else{
					removeDuplicates = true;
					cerr<<"Warning: "<<item<<" is not a valid option for RemoveDuplicates; it has to be TRUE or FALSE"<<endl;
					cerr<<"Warning: RemoveDuplicates is set to TRUE"<<endl;
				}
			}
            else if (item.compare("Type")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("chip-seq")==0 || tempStr.compare("chipseq")==0 || tempStr.compare("chip")==0)
					isCHIPseq = true;
				else{
					isCHIPseq = false;
					cout<<"Considering that the experiment is DNAse or ATAC-seq"<<endl;
				}
			}
			else if (item.compare("PrintPosterior")==0){
	                        transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
	                        if (tempStr.compare("true")==0)
	                                posterior = true;
	                        else if (tempStr.compare("false")==0)
	                                posterior = false;
	                        else{
	                                posterior = false;
	                                cerr<<"Warning: "<<item<<" is not a valid option for printPosterior; it has to be TRUE or FALSE"<<endl;
	                                cerr<<"Warning: printPosterior is set to FALSE"<<endl;
	                        }
	                }
            else if (item.compare("pairedEnds")==0){
	                        transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
	                        if (tempStr.compare("true")==0)
	                                pairedEnds = true;
	                        else if (tempStr.compare("false")==0)
	                                pairedEnds = false;
	                        else{
	                                pairedEnds = false;
	                                cerr<<"Warning: "<<item<<" is not a valid option for pairedEnds; it has to be TRUE or FALSE"<<endl;
	                                cerr<<"Warning: pairedEnds is set to FALSE"<<endl;
	                        }
            }
			else if (item.compare("PosteriorProb")==0)
				posterior_threshold = atof(tempStr.c_str());
			else if(item.compare("blackListFile")== 0)
				blacklist = tempStr;

			else if (item.compare("PrintBedGraph")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("true")==0)
					bedgraph = true;
				else if (tempStr.compare("false")==0)
					bedgraph = false;
				else{
					bedgraph = false;
					cerr<<"Warning: "<<item<<" is not a valid option for printBedGraph; it has to be TRUE or FALSE"<<endl;
					cerr<<"Warning: "<<item<<"is set to FALSE"<<endl;
				}
			}

			else if (item.compare("CallPeaks")==0){
				transform(tempStr.begin(), tempStr.end(),tempStr.begin(), ::tolower);
				if (tempStr.compare("true")==0)
					callpeaks = true;
				else if (tempStr.compare("false")==0)
					callpeaks = false;
				else{
					callpeaks = true;
					cerr<<"Warning: "<<item<<" is not a valid option for CallPeaks; it has to be TRUE or FALSE"<<endl;
					cerr<<"Warning: "<<item<<"is set to True"<<endl;
				}
			}

			else
				cerr<<"Warning: "<<item<<" is not a parameter for HMCan"<<endl;
		}

		// validating parameters
		if (GC_index=="")
		{
			cerr<<"Please provide a valid file name for GC content information file"<<endl;
			exit(EXIT_FAILURE);
		}

		if (path=="")
		{
			cerr<<"Please provide a valid path for your target genome"<<endl;
			exit(EXIT_FAILURE);
		}
        if (pairedEnds) {
            if (t == BED) {
                cerr<<"Warning: cannot use Paired End option on BED input files"<<endl;
                pairedEnds=false;
            } else {
                min_length = 0;
                median_length = 0;
                max_length = 0;
                cout << "Paired End mode is enabled. It will not take into account provided values for fragment lengths; will evaluate fragment sizes from the provided SAM/BAM files\n";
            }
        }
        if (!pairedEnds) {
            if(min_length<1){
                min_length = 145;
                cerr<<"Warning: Invalid parameter, minimum fragment length is set to its default value 145"<<endl;
            }

            if(median_length<1){
                median_length = 150;
                cerr<<"Warning: Invalid parameter, median fragment length is set to its default value 150"<<endl;
            }

            if(max_length<1){
                max_length = 155;
                cerr<<"Warning: Invalid or missed parameter, maximum fragment length is set to its default value 155"<<endl;
            }
        }
        if (!isCHIPseq) {
            max_length = 70;
            median_length = 50;
            min_length = 20;
            cout<<"Warning: Overwriting parameters of read extension to 20/50/60"<<endl;
        }

		if(bin_size<1){
			bin_size = 50;
			cerr<<"Warning: Invalid or missed parameter, bin length length is set to its default value 50"<<endl;
		}

		if(large_bin_size<1){
			large_bin_size = 100000;
			cerr<<"Warning: Invalid or missed parameter, large bin size length is set to its default value 100000"<<endl;
		}

		if(pvalue_threshold<0 || pvalue_threshold>1){
			pvalue_threshold = 0.01;
			cerr<<"Warning: Invalid or missed parameter, P-value threshold is set to its default value 0.01"<<endl;
		}

		if(merge_dist<0){
			 merge_dist= 200;
			cerr<<"Warning: Invalid or missed parameter, merge distance is set to its default value 200"<<endl;
		}

		if(iter_threshold<0){
			iter_threshold = 5;
			cerr<<"Warning: Invalid or missed parameter, iterations score threshold is set to its default value 5"<<endl;
		}

		if(final_threshold<0){
			final_threshold = 0;
			cerr<<"Warning: Invalid or missed parameter, final score threshold is set to its default value 0"<<endl;
		}

		if (max_iter<1){
			max_iter = 10;
			cerr<<"Warning: Invalid or missed parameter, maximum iterations is set to its default value 10"<<endl;
		}
		if(posterior_threshold == -1){
			posterior_threshold = 0.7;
			cerr<<"Warning: Invalid or missed parameter, PosteriorProb is set to its default value 0.7"<<endl;
		}



}


void Parser::print(){
	cout<<GC_index<<endl;
	cout<<path<<endl;
	cout<<min_length<<endl;
	cout<<median_length<<endl;
	cout<<max_length<<endl;
	cout<<large_bin_size<<endl;
	cout<<bin_size<<endl;
	cout<<pvalue_threshold<<endl;
	cout<<merge_dist<<endl;
	cout<<iter_threshold <<endl;
	cout<<final_threshold <<endl;
	cout<<max_iter<<endl;
	cout<<blacklist<<endl;
	cout << pairedEnds<<endl;
	cout << removeDuplicates<<endl;
    cout << GC_merge_dist<<endl;
    cout << calculateEmpiricalPvalue<<endl;
    cout << ifKeepChIPdensityIntact<<endl;
    cout << CNAnormalization<<endl;

}
