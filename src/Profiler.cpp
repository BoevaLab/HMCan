/*************************************************************************
Copyright (c) 2013, Haitham ASHOOR

>>> SOURCE LICENSE >>>
r
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
#include "Profiler.h"
#include "linreg.h"
#include "segmentation.h"
#include "freec.h"
#include "poissondistr.h"
#include "utils.h"

#include <sstream>
#include <set>
using namespace std;


Profiler::Profiler(const char * GC_index, int left, int med, int right,int bin_length,int large_bin_size,int merge_dist,float pvalue_threshold,
                bool pairedEnds, bool isChIPseq, bool hasControl, int GC_merge_dist, bool calculateEmpiricalPvalue) {

	ifstream gc;
	string chr;
	gc.open(GC_index);

	if(!gc){
		cerr<<"Error:can not open GC information file"<<endl;
		exit(1);
	}
	read_gc_profile(gc);
	gc.close();
	this->right = right;
	this->left = left;
	this->med = med;
	this->bin_length = bin_length;
	this->large_bin_size = large_bin_size;
	this->merge_dist = merge_dist;
	this->pvalue_threshold = pvalue_threshold;
	pairedEnds_=pairedEnds;
	calculate_density_coefs();
	min_obs = 0;
	min_obs_Input=0;
	max_obs = 0;
    max_obs_Input=0;
	isChIPseq_=isChIPseq;
    hasControl_=hasControl;
    GC_mergeDist=GC_merge_dist;
    calculateEmpiricalPvalue_=calculateEmpiricalPvalue;
}




Profiler::~Profiler() {

	map<string,float*>::iterator chr_it;
	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.begin();++chr_it)
		delete (*chr_it).second;

	for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.begin();++chr_it)
			delete (*chr_it).second;

	for(chr_it=sampled_GC.begin();chr_it!=sampled_GC.begin();++chr_it)
			delete (*chr_it).second;

     for(chr_it=peaks.begin();chr_it!=peaks.end();++chr_it)
    	 delete (*chr_it).second;
}

void Profiler::build_single_profile(vector<DNA_fragment >& tags, string chr, bool targetORcontrol, bool CNAnormalization){

	// allocate memory for density
	if (targetORcontrol)
		sampled_target_density[chr] =  new float[sampled_sizes[chr]];
	else
		sampled_control_density[chr] =  new float[sampled_sizes[chr]];

    if (targetORcontrol==0 && calculateEmpiricalPvalue_) { //need to calculate empirical p-value
            sampled_control_density1[chr]=  new float[sampled_sizes[chr]];
            sampled_control_density2[chr]=  new float[sampled_sizes[chr]];
    }


	for(int i=0;i<sampled_sizes[chr];i++) {
		if (targetORcontrol)
			sampled_target_density[chr][i] = 0;
		else
			sampled_control_density[chr][i] = 0;
        if (targetORcontrol==0 && calculateEmpiricalPvalue_) {
			sampled_control_density1[chr][i] = 0;
			sampled_control_density2[chr][i] = 0;
        }
	}

	for(vector<DNA_fragment>::iterator fragment_it=tags.begin(); fragment_it !=tags.end();++fragment_it ){
        int flip=0;
        if (calculateEmpiricalPvalue_ && targetORcontrol==0) { // flip a coin
                flip = (rand()%2) ;
        }
        if ((*fragment_it).fragmentLength==0 || (*fragment_it).fragmentLength>right && isChIPseq_) { //single end or should be considered as SE due to a strange mapping:
            int fragment_start,fragment_left,fragment_med, fragment_right;
            if ((*fragment_it).strand==0) { //forward read:
                fragment_start = (*fragment_it).start;
                fragment_left  = (*fragment_it).start+left-1>sizes[chr] ? sizes[chr]-1:(*fragment_it).start+left-1;
                fragment_med = 	(*fragment_it).start+med-1>sizes[chr]?sizes[chr]-1:(*fragment_it).start+med-1;
                fragment_right = (*fragment_it).start+right-1>sizes[chr]? sizes[chr]-1:(*fragment_it).start+right-1;


                //start at multiple of bins position
                int start_point = fragment_start;
                while ((start_point+1)%bin_length!=0)
                    start_point++;

                //calculate density using triangular distribution
                //check the FindPeaks program for more info
                for(int i=start_point;i<=fragment_right;i+=bin_length){
                    int index = ((i+1)/bin_length)-1;
                    if(index>=sampled_sizes[chr]){
                        cerr<<"Incorrect reference genome, please check your reference genome"<<endl;
                        exit(0);
                    }
                    int x = i-fragment_start;
                    float point_density;
                    if (i<=fragment_left) {
                        point_density=1;
                    } else if (i<=fragment_med)   {
                        float hx = slope1 * x + b1;
                        point_density= 1 - ((hx * (x-left)) /2);
                    } else {
                        float hx = slope2 * x + b2;
                        point_density = (1 - left_area)-(right_area - (hx * (float)(right-x) /2));
                    }
                    if(targetORcontrol)
                        sampled_target_density[chr][index]+=point_density;
                    else
                        sampled_control_density[chr][index]+=point_density;

                    if(targetORcontrol==0 && calculateEmpiricalPvalue_ && flip) {
                        sampled_control_density1[chr][index]+=point_density;
                    }
                    if(targetORcontrol==0 && calculateEmpiricalPvalue_ && !flip) {
                        sampled_control_density2[chr][index]+=point_density;
                    }
                }
            } else { //reverse read:
                fragment_start = (*fragment_it).start+right-1>sizes[chr]? sizes[chr]-1:(*fragment_it).start+right-1; //rightmost point
                fragment_left  = (*fragment_it).start; //leftmost point
                fragment_med = 	fragment_start - med;
                fragment_right = fragment_start - left;
                if (fragment_med<0)
                    fragment_med=0;
                if (fragment_right<0)
                    fragment_right=0;

                 //start at multiple of bins position
                int start_point = fragment_left;
                while ((start_point+1)%bin_length!=0)
                    start_point++;

                //calculate density using triangular distribution
                //check the FindPeaks program for more info
                for(int i=start_point;i<=fragment_start;i+=bin_length){
                    int index = ((i+1)/bin_length)-1;
                    if(index>=sampled_sizes[chr]){
                        cerr<<"Incorrect reference genome, please check your reference genome"<<endl;
                        exit(0);
                    }
                    int x = fragment_start-i;
                    float point_density;
                    if (i>=fragment_right) {
                        point_density=1;
                    } else if (i>=fragment_med)   {
                        float hx = slope1 * x + b1;
                        point_density = 1 - ((hx * (x-left)) /2);

                    } else {
                        float hx = slope2 * x + b2;
                        point_density = (1 - left_area)-(right_area - (hx * (float)(right-x) /2));
                    }
                    if(targetORcontrol)
                        sampled_target_density[chr][index]+=point_density;
                    else
                        sampled_control_density[chr][index]+=point_density;

                    if(targetORcontrol==0 && calculateEmpiricalPvalue_ && flip) {
                        sampled_control_density1[chr][index]+=point_density;
                    }
                    if(targetORcontrol==0 && calculateEmpiricalPvalue_ && !flip) {
                        sampled_control_density2[chr][index]+=point_density;
                    }
                }
            }
        } else { // fragmentLength >0, PE:
            int fragment_start, fragment_right;

            fragment_start = (*fragment_it).start;
            fragment_right = (*fragment_it).start+(*fragment_it).fragmentLength-1>sizes[chr]? sizes[chr]-1:(*fragment_it).start+(*fragment_it).fragmentLength-1;

            //start at multiple of bins position
            int start_point = fragment_start;
            while ((start_point+1)%bin_length!=0)
                start_point++;

            //add 1 to density
            for(int i=start_point;i<=fragment_right;i+=bin_length){
                int index = ((i+1)/bin_length)-1;
                if(index>=sampled_sizes[chr]){
                    cerr<<"Incorrect reference genome, please check your reference genome"<<endl;
                    exit(0);
                }
                if(targetORcontrol)
                    sampled_target_density[chr][index]+= 1;
                else
                    sampled_control_density[chr][index]+=1;
                if(targetORcontrol==0 && calculateEmpiricalPvalue_ && flip) {
                        sampled_control_density1[chr][index]+=1;
                }
                if(targetORcontrol==0 && calculateEmpiricalPvalue_ && !flip) {
                        sampled_control_density2[chr][index]+=1;
                }

            }
        }
	}

	if (pairedEnds_ && isChIPseq_) { //divide density by 2 for PE datasets:
        for(int i=0;i<sampled_sizes[chr];i++){


            if(targetORcontrol) {
                sampled_target_density[chr][i] =sampled_target_density[chr][i]/2.0;
            } else {
                sampled_control_density[chr][i]=sampled_control_density[chr][i]/2.0;
            }
            if(targetORcontrol==0 && calculateEmpiricalPvalue_) {
                sampled_control_density1[chr][i]=sampled_control_density1[chr][i]/2.0;
                sampled_control_density2[chr][i]=sampled_control_density2[chr][i]/2.0;
            }
        }
	}

	if (CNAnormalization) {
        //correct for copy number
        for(int i=0;i<sampled_sizes[chr];i++){
            int index = i*bin_length+bin_length;
            int bin_index = index/large_bin_size;
            float median = medians[chr][bin_index];
            if(median>0.3){ // check for low copy number if found assign unknown label for that
                if (targetORcontrol){
                    sampled_target_density[chr][i] =sampled_target_density[chr][i]/median;
                }else {
                    sampled_control_density[chr][i]=sampled_control_density[chr][i]/median;
                }
                if(targetORcontrol==0 && calculateEmpiricalPvalue_ ) {
                    sampled_control_density1[chr][i]=sampled_control_density1[chr][i]/median;
                    sampled_control_density2[chr][i]=sampled_control_density2[chr][i]/median;
                }
            }

        }
    }

}

void Profiler::build_profile_wig(const char * genomeLengthFile, const char * wigfilename){

    ifstream genomeLengthFileStream(genomeLengthFile);
    string line;
    std::vector<std::string> strs;
    map<string,float*>::iterator chr_it;

    if (!genomeLengthFileStream.is_open()) {
        cerr << "ERROR: Unable to open " << genomeLengthFile << " file for reading\n";
        cerr << "Please provide in the config file .fa.fai file with the reference genome chromosome lengths (e.g. hg19.fa.fai)\n";
        cerr << "Use option \"genomeLengthFile yourPath/hg19.fa.fai\"\n";
        exit (-1);
    }
    else {
        while (!genomeLengthFileStream.eof()) {
            getline(genomeLengthFileStream, line);
            if (line.length()>0) {
                strs = split(line, '\t');
                sampled_sizes[strs[0].c_str()] = atoi(strs[1].c_str())/bin_length+1; //bin_length is the span size in wig file
            }
        }
    }
    strs.clear();
    genomeLengthFileStream.close();

    // Read the wig file.
    ifstream infile;
    infile.open(wigfilename);
    if (!infile){
        cerr<<"Error: "<<wigfilename<<" can not be found"<<endl;
        exit(1);
    }

    std::string chr;
    std::string chrn;
    int ind;
    chr = "chr0";
    bool isFixedStep = true;

	while(!infile.eof()){
        getline(infile,line);
        if(line.compare("")==0)
            continue;
        istringstream stream(line);
        vector<string> tokens;
        copy(istream_iterator<string>(stream),
			 istream_iterator<string>(),
			 back_inserter<vector<string> >(tokens)); // split the line
        if (tokens[0].compare("track")==0)
            continue;   //skip the track heade
        else if (tokens[0].compare("variableStep")==0) {
            isFixedStep = false;
            strs = split(tokens[1], '=');
            chrn = strs[1].c_str();
            if (chrn.compare(chr)==0)
                continue;
            else {
                chr = strs[1].c_str();
                //cout<<"Building Density profile for chromosome "<<chr<<"......."<<endl;
                sampled_target_density[chr] =  new float[sampled_sizes[chr]];
                for(int i=0;i<sampled_sizes[chr];i++)
                    sampled_target_density[chr][i] = 0;
            }
        } else if (tokens[0].compare("fixedStep")==0) {
            strs = split(tokens[1], '=');
            chrn = strs[1].c_str();
            int startPosition = 1;
            if (tokens.size() > 2 && (split(tokens[2], '=')[0]).compare("start")==0)
                startPosition = atoi((split(tokens[2], '=')[1]).c_str());
            else if (tokens.size() > 3 && (split(tokens[3], '=')[0]).compare("start")==0)
                startPosition = atoi((split(tokens[3], '=')[1]).c_str());

            ind = startPosition/bin_length;
            if ( (tokens.size() > 2 && (split(tokens[2], '=')[0]).compare("step")==0 && atoi((split(tokens[2], '=')[1]).c_str())!=bin_length)
             || (tokens.size() <=2 && bin_length!=1)
             || (tokens.size() > 3 && (split(tokens[3], '=')[0]).compare("step")==0 && atoi((split(tokens[3], '=')[1]).c_str())!=bin_length)) {
                cerr << "WARNING: the bin size of your wig file is not equal to the small bin length provided in the config file..\n";
                cerr << "Please change the value of smallBinLength in the config file to the right value of step (from the .wig file)" <<endl;
                cerr << "Exiting..";
                exit(0);
            }

            if (chrn.compare(chr)==0)
                continue;
            else {
                chr = strs[1].c_str();
                //cout<<"Building Density profile for chromosome "<<chr<<"......."<<endl;
                sampled_target_density[chr] =  new float[sampled_sizes[chr]];
                for(int i=0;i<sampled_sizes[chr];i++)
                    sampled_target_density[chr][i] = 0;
            }
        }
        else {
            if (!isFixedStep) {
                ind = (atof(tokens[0].c_str())/bin_length); // Elnaz, you do not need "+1", I think;
                sampled_target_density[chr][ind] = atof(tokens[1].c_str());
            } else {
                sampled_target_density[chr][ind] = atof(tokens[0].c_str());
                ind++;
            }
        }
    }

    for (chr_it = sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		peaks[(*chr_it).first] = new float[sampled_sizes[(*chr_it).first]];
	}

	call_peaks(GC_mergeDist);
	cout<<"call_peaks done...."<<endl;
	derive_transition_probabilities();
	cout<<"derive_transition_probabilities done...."<<endl;
	generate_observations();
	cout<<"generate_observations done...."<<endl;
	count_emissions();
	cout<<"count_emissions done...."<<endl;

}


void Profiler::build_profile(map<string,vector<DNA_fragment > >& data, map<string,vector<DNA_fragment > >& control, string& path,
		vector<BedEntry>& blacklist, bool keepChIPdensityIntact, bool CNAnormalization){

	float N=0,M=0;
	map<string,vector<DNA_fragment > >::iterator chr_it;


	clean_reads(data,control);
	for(chr_it = data.begin();chr_it!=data.end();++chr_it){
		N+=(*chr_it).second.size();
	}


	for(chr_it = control.begin();chr_it!=control.end();++chr_it)
		M+=(*chr_it).second.size();

    if (M!=0 && N!=0) {
        reads_ratio = M/N;
    } else {
        reads_ratio=1;
    }

	cout<<"Reads ratio is: "<<reads_ratio<<endl;
	for(chr_it = data.begin();chr_it!=data.end();++chr_it){
		if ((*chr_it).first.compare("chrM")==0)
			continue;
		string chromosome_path;
		if (path[path.size()-1] == '/')
			chromosome_path = path+(*chr_it).first+".fa";
		else
			chromosome_path = path+"/"+(*chr_it).first+".fa";

		if (right==0 && isChIPseq_) { //PE data, first chromosome being read:
            //evaluate maximum, minimum and median fragment length:
            calculateFragLengthDist((*chr_it).second);
            calculate_density_coefs();

		}
		string chr_seq = Read_chr(chromosome_path.c_str());
		sizes[(*chr_it).first] = chr_seq.size();
		cout<<(*chr_it).first<<"\t"<<chr_seq.size()<<endl;
		sampled_sizes[(*chr_it).first]=sizes[(*chr_it).first]/bin_length+1;
		extend_reads((*chr_it).second,sizes[(*chr_it).first]);
		extend_reads(control[(*chr_it).first],sizes[(*chr_it).first]);

		calculate_GC_sampled_bins((*chr_it).first,chr_seq);
	}

	if (CNAnormalization) {
        cout<<"calculating copy number alternations........"<<endl;
        if (hasControl_)
            call_freec(control, true); //use Input
        else {
            call_freec(data, false); //use ChIP or ATAC or DNAse, so will need to remove signal regions first to get a smooth signal
        }
    } else {
        cout<<"HMCan is not going to correct for copy number alternations........"<<endl;
    }

	for (chr_it = data.begin();chr_it!=data.end();++chr_it){
		cout<<"Building Density profile for chromosome "<<(*chr_it).first<<"......."<<endl;
		build_single_profile((*chr_it).second,(*chr_it).first,true,CNAnormalization);
		if (hasControl_) build_single_profile(control[(*chr_it).first],(*chr_it).first,false,CNAnormalization);
	}

	remove_blacklist(blacklist,sampled_target_density);
	if (hasControl_) remove_blacklist(blacklist,sampled_control_density);
	if (hasControl_ && calculateEmpiricalPvalue_) {
        remove_blacklist(blacklist,sampled_control_density1);
        remove_blacklist(blacklist,sampled_control_density1);
	}

	for (chr_it = data.begin();chr_it!=data.end();++chr_it){
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++)
			if (sampled_target_density[(*chr_it).first][i]!=-1){
				if (reads_ratio>1 && !keepChIPdensityIntact)
					sampled_target_density[(*chr_it).first][i] = sampled_target_density[(*chr_it).first][i]*reads_ratio;
				else
					if (hasControl_) sampled_control_density[(*chr_it).first][i] = sampled_control_density[(*chr_it).first][i]/reads_ratio;
			}
		peaks[(*chr_it).first] = new float[sampled_sizes[(*chr_it).first]];
		if (calculateEmpiricalPvalue_)
            peaks_Input[(*chr_it).first] = new float[sampled_sizes[(*chr_it).first]];

	}
    if (isChIPseq_){
        cout<<"Normalize for GC content bias and noise ratio........."<<endl;
        call_peaks(GC_mergeDist); // XXX
        Normalize_GC();
    } else {
        cout<<"Will not normalize for GC-content bias and noise ratio........."<<endl;
    }
	call_peaks(merge_dist);
	derive_transition_probabilities();
	generate_observations();
	count_emissions();

	if (calculateEmpiricalPvalue_) {
        call_peaks_Input(merge_dist);
        derive_transition_probabilities_Input();
        generate_observations_Input();
        count_emissions_Input();
	}

}

int Profiler::getMaxLength() {
    return right;
}

void Profiler::read_gc_profile(ifstream& gc_profile){
	string line;
	while(getline(gc_profile,line)){
		vector<string> tokens;
		istringstream stream(line);
		copy(istream_iterator<string>(stream),
		     istream_iterator<string>(),
			 back_inserter<vector<string> >(tokens));
		if (tokens.size()<5){
			cerr<<"Warning: the line: "<<line<<" ,does not follow format. Line will be ignored"<<endl;
                        continue;
		}
		string key = tokens[0];
		float content = atof(tokens[2].c_str());
		float mapability = atof(tokens[4].c_str());
		GC_profile[key].push_back(content);
		notNprofile[key].push_back(mapability);

		}
}


void Profiler::call_freec(map<string,vector<DNA_fragment > >& tags, bool isInput){
	map <string, vector <float> > read_count;
	//map <string, vector <float> > ratio_profile;
	map <string, vector <int> > breakpoints;
	int read_counts = 0;
    float breakpointsointThreshold = 0.8;
	if (isInput) {
        for (map<string,vector<DNA_fragment > >::iterator chr_it = tags.begin();chr_it!=tags.end();++chr_it) {
            read_counts+=(*chr_it).second.size();
			int bins_count =GC_profile[(*chr_it).first].size();
			read_count[(*chr_it).first] = vector<float>(bins_count,0);
			for(vector<DNA_fragment>::iterator v = (*chr_it).second.begin();v!=(*chr_it).second.end();++v) {
				int bin_index = (*v).start/large_bin_size;
				if (bin_index >= bins_count){
					cerr<<"Error:" <<"There is a problem with "<<(*chr_it).first<<" size please check your GC Index file!"<<endl;
                    cerr<<"Check that 'largeBinLength' corresponds to the bin size in 'GCIndex'"<<endl;
					exit(1);
				}
				read_count[(*chr_it).first][bin_index]+=1;
			}
		}
	} else {
        breakpointsointThreshold=1.6; //as the data will be much noisier
        for (map<string,vector<DNA_fragment > >::iterator chr_it = tags.begin();chr_it!=tags.end();++chr_it) {
            vector<DNA_fragment> cleanTags ( (*chr_it).second);
            read_counts+=cleanTags.size();
            int localBinSize = 1000;
            if (large_bin_size<localBinSize) {
                localBinSize = large_bin_size/4; // should not happen
                cout << "Warning: please increase the value of largeBinLength to at least 4 Kb\n";
            } else {
                //check that large_bin_size is multiple of localBinSize
                if (large_bin_size % localBinSize != 0) {
                    localBinSize=500;
                    cout << "Warning: please set the value of largeBinLength multiple of 1000 bp\n";
                }
            }
            //now for each large_bin find local bins to remove:
            float quantileToRemove = 1-ATACSEQPEAKPROP;
            if (isChIPseq_) quantileToRemove=1-CHIPSEQPEAKPROP;
            cout <<"Estimating copy numbers without a control sample. Will remove "<<(100-quantileToRemove*100)<< "% of windows that can potentially contain an ATAC/ChIP/DNAse signal\n";
            int numberOfsmallBinsInBig=large_bin_size/localBinSize;
            int bins_count =GC_profile[(*chr_it).first].size();
            int localBinCount =numberOfsmallBinsInBig*bins_count +1;
			vector<float> localCounts = vector<float>(localBinCount,0);
			for(vector<DNA_fragment>::iterator v = cleanTags.begin();v!=cleanTags.end();++v) {
				int bin_index = (*v).start/localBinSize;
				if (bin_index >= localBinCount){
					cerr<<"Error:" <<"There is a problem with "<<(*chr_it).first<<" size please check your GC Index file!"<<endl; // should never get here if parameters are correct
                    cerr<<"Check that 'largeBinLength' corresponds to the bin size in 'GCIndex'"<<endl;
					exit(1);
				}
				localCounts[bin_index]+=1;
			}
            read_count[(*chr_it).first] = vector<float>(bins_count,0);
            for (int bin_index=0; bin_index<bins_count; bin_index++) {
                vector<float> localCountsInThisBin(localCounts.begin() + bin_index*numberOfsmallBinsInBig,localCounts.begin() + bin_index*numberOfsmallBinsInBig+numberOfsmallBinsInBig);
                sort( localCountsInThisBin.begin(), localCountsInThisBin.end());
                int sum_of_elems = 0; int winCount=0;
                int maxElementToTake=floor(localCountsInThisBin.size()*quantileToRemove);

//                float IQRby = localCountsInThisBin[floor(localCountsInThisBin.size()*0.75)]-localCountsInThisBin[floor(localCountsInThisBin.size()*0.25)];
//                float maxElementValueToKeep = IQRby+localCountsInThisBin[floor(localCountsInThisBin.size()*0.5)];

                vector<float>::iterator nth = localCountsInThisBin.begin() + maxElementToTake;
                for(std::vector<float>::iterator it = localCountsInThisBin.begin(); it !=nth; ++it) {
                    sum_of_elems += *it;winCount++;
                }
//                vector<float>::iterator nth = localCountsInThisBin.begin() + maxElementToTake;
//                for(std::vector<float>::iterator it = localCountsInThisBin.begin(); *it <= maxElementValueToKeep; ++it) {
//                    sum_of_elems += *it;winCount++;
//                }

                if(winCount>0)
                    read_count[(*chr_it).first][bin_index]=sum_of_elems; // or sum_of_elems*numberOfsmallBinsInBig/winCount; if uses maxElementValueToKeep
                localCountsInThisBin.clear();
            }
            localCounts.clear();
			cleanTags.clear();

            //clean further the copy number profiles:
            removeOutliers(read_count[(*chr_it).first],1.4,5);
		}
	}
	cout<<"Total reads count read by FREEC is: "<<read_counts<<endl;
	int ploidy = 4; // to cover all cases; should not affect negatively cases with ploidy==2 or 3
	recalculateRatioUsingCG ( ploidy, read_count, GC_profile, notNprofile,ratio_profile);
	calculateBreakpoints(ratio_profile, breakpoints, breakpointsointThreshold);
	calculateCopyNumberMedians(ploidy,breakpoints,medians,ratio_profile);
}


void Profiler::calculateFragLengthDist(std::vector<DNA_fragment >& tags) {
    vector<int> fragmentLengths;
    for(vector<DNA_fragment>::iterator tag_it = tags.begin(); tag_it!=tags.end();++tag_it)
		if ((*tag_it).fragmentLength >0) {
            fragmentLengths.push_back((*tag_it).fragmentLength);
        }
    if (fragmentLengths.size()>0) {
        cout << "Fragment length distribution will be evaluated based on "<<fragmentLengths.size()<<" fragments\n";
    }else {
        cerr << "Error: Could not find paired-end reads to evaluate fragment size\n";
        cerr << "Remove 'pairedEnds True' option from your config file and provide fragment length distribution\n";
        exit(1);
    }
    getFragmentLengthDistr(fragmentLengths);
    cout << "Evaluated fragment sizes: "<<left<<", " << med << ", "<< right <<endl;
    return;
}

void Profiler::getFragmentLengthDistr (std::vector<int> & myvector) {
    vector<int> data (myvector);
    int ndatapoints = data.size();
    if (ndatapoints==0) {
        cerr << "Error: zero values to calculate medians..\n";
        exit(-1);
    }
    sort(data.begin(),data.end());
    // Get median
    med = ndatapoints % 2 == 1 ? data[(ndatapoints-1)/2] : (data[ndatapoints/2 - 1] + data[ndatapoints/2])/float(2.0);
    int upperboundary=med*10;
    int maxIndexToKeep = 0;
    for (; maxIndexToKeep<ndatapoints;maxIndexToKeep++) {
        if(data[maxIndexToKeep]>upperboundary) {
            break;
        }
    }
    data.resize(maxIndexToKeep);
    right=data[ceil(float((maxIndexToKeep-1)*99.0)/100)];
    left=data[ceil(float((maxIndexToKeep-1)*5.0)/100)];

    data.clear();

    return ;
}

void Profiler:: extend_reads(std::vector<DNA_fragment >& tags, int chr_size){

    if (isChIPseq_) {
        for(vector<DNA_fragment>::iterator tag_it = tags.begin(); tag_it!=tags.end();++tag_it) {
            if ((*tag_it).strand == 1) {
                if ((*tag_it).fragmentLength>0) {
                    if ((*tag_it).start-(*tag_it).fragmentLength-1<0)
                        (*tag_it).start = 0;
                    else
                        (*tag_it).start = (*tag_it).start-(*tag_it).fragmentLength+1;
                } else {
                    if ((*tag_it).start-right-1<0)
                        (*tag_it).start = 0;
                    else
                        (*tag_it).start = (*tag_it).start-right+1;
                }
            }
        }
    } else{
        int nucleLenThreshold=130;
        vector<DNA_fragment> newTags;
        //int right = 70; should be already set
        for(vector<DNA_fragment>::iterator tag_it = tags.begin(); tag_it!=tags.end();++tag_it) {
            if ((*tag_it).strand == 1) { //reverse read
                if ((*tag_it).fragmentLength >= nucleLenThreshold) {
                    (*tag_it).fragmentLength=0; //cut the fragment
                }
                //create a mirrow fragment:
                DNA_fragment read;read.fragmentLength=0; read.strand=0; read.start=(*tag_it).start+1; newTags.push_back(read);
                if ((*tag_it).fragmentLength>0) {
                    if ((*tag_it).start-(*tag_it).fragmentLength-1<0)
                        (*tag_it).start = 0;
                    else
                        (*tag_it).start = (*tag_it).start-(*tag_it).fragmentLength+1;
                    //and add a left tail read:
                    DNA_fragment read;read.fragmentLength=0; read.strand=1; read.start=(*tag_it).start-right; if (read.start<0)read.start=0; newTags.push_back(read);
                } else {
                    if ((*tag_it).start-right-1<0)
                        (*tag_it).start = 0;
                    else
                        (*tag_it).start = (*tag_it).start-right+1;
                }
            } else {  //forward read
                if ((*tag_it).fragmentLength >= nucleLenThreshold) {
                    (*tag_it).fragmentLength=0; //cut the fragment
                }
                //create a mirrow fragment:
                DNA_fragment read;read.fragmentLength=0; read.strand=1; read.start=(*tag_it).start-right; if (read.start<0)read.start=0; newTags.push_back(read);
                if ((*tag_it).fragmentLength >0) {
                    //create a reads on the right of the fragment:
                    DNA_fragment read;read.fragmentLength=0; read.strand=0; read.start=(*tag_it).start+(*tag_it).fragmentLength+1; newTags.push_back(read);
                }
            }
        }
        //add new tags:
        tags.insert( tags.end(), newTags.begin(), newTags.end() );
    }



	//remove out of range items
	vector<DNA_fragment>::iterator it = tags.begin();

	while(it!=tags.end()){
		if ((*it).start >= chr_size){
			it = tags.erase(it);
		}
		else{
			++it;
		}
	}

    sort(tags.begin(),tags.end(),compare); // resort to take extension into consideration

}


void Profiler::calculate_density_coefs(){
	h = 2/(float)(right-left);
	slope1 = h/(float)(med-left);
	slope2 = -h/(float)(right-med);
	b1 = -slope1*left;
	b2 = -slope2*right;
	left_area = (h*(float)(med-left))/2;
	right_area = (h*(float)(right-med)/2);

}




void Profiler::Normalize_GC(){

	float chip_accum_density[101],chip_windows_count[101],control_accum_density[101],control_windows_count[101],
			chip_percent[101],control_percent[101],GC_stratum=0,chip_lambda=0,control_lambda=0,max_chip=0,
			max_control=0;
	int startum_windows;
	vector<vector<float> > density;
	map<string,float*>::iterator chr_it;
	vector<float> chip,control;


	//normalize target
	for(int i=0;i<101;i++){
		chip_accum_density[i]=0;
		chip_windows_count[i]=0;
		density.push_back(vector<float>());
	}


	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++)
			if (peaks[(*chr_it).first][i]==0&& (*chr_it).second[i]!=-1 && sampled_GC[(*chr_it).first][i]!=-1){
				int index = iround(sampled_GC[(*chr_it).first][i]); //(int)(sampled_GC[(*chr_it).first][i]*100);
				chip_windows_count[index]++;
				density[index].push_back((*chr_it).second[i]);
			}
	}

	for(int i=0;i<101;i++)
        if(chip_windows_count[i]>0){
				chip_windows_count[i]=0;
				int percent_to_remove = density[i].size()/10; //10%
				if (!isChIPseq_) {
                    percent_to_remove = density[i].size()/6; //15%
				}
				sort(density[i].begin(),density[i].end());
				float sum=0;

  				for(unsigned int j=percent_to_remove;j<density[i].size()-percent_to_remove;j++){
					sum+=density[i][j];
					chip_windows_count[i]++;
				}
				chip_accum_density[i] = sum;
				if (chip_windows_count[i]<10) {chip_windows_count[i]=0;chip_accum_density[i] =0;} //ignore GC-bin that are very rare
        }

//    startum_windows=0;
    int minGCPercentoConsider=27; //should be below 35 but >0
    int maxGCPercentoConsider=75; //should be above 66 but <100

//    for(int i=0;i<27;i++){
//        GC_stratum+=chip_accum_density[i];
//        startum_windows+=chip_windows_count[i];
//    }


    for(int i=(minGCPercentoConsider+1);i<=35;i+=2){
			GC_stratum=(chip_accum_density[i]+chip_accum_density[i-1])/(chip_windows_count[i]+chip_windows_count[i-1]);
			chip_percent[i] = chip_percent[i-1] = GC_stratum;
    }

    for(int i=0;i<minGCPercentoConsider;i++)
        chip_percent[i] =chip_percent[minGCPercentoConsider]; // used to be until v1.43: =GC_stratum/startum_windows;


    for(int i=35;i<66;i++){
        if (chip_windows_count[i]!=0)
            chip_percent[i]=chip_accum_density[i]/chip_windows_count[i];
        else {
            chip_percent[i]=chip_percent[i-1];
        }
    }



    for(int i = 66;i<maxGCPercentoConsider;i+=2){
        GC_stratum=(chip_accum_density[i]+chip_accum_density[i-1])/(chip_windows_count[i]+chip_windows_count[i-1]);
        chip_percent[i] = chip_percent[i-1] = GC_stratum;
    }

//    GC_stratum=0;
//    startum_windows=0;
//    for(int i=maxGCPercentoConsider;i<101;i++){
//			GC_stratum+=chip_accum_density[i];
//			startum_windows+=chip_windows_count[i];
//    }

    for(int i=maxGCPercentoConsider;i<101;i++)
        chip_percent[i] = chip_percent[maxGCPercentoConsider-1];//used to be until v1.43: =GC_stratum/startum_windows;

	density.clear();



								//normalize for control data
	//*****************************************************************************************************************************
 if (hasControl_) {

	for(int i=0;i<101;i++){
		control_accum_density[i]=0;
		control_windows_count[i]=0;
		density.push_back(vector<float>());
	}

	for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.end();++chr_it){
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++)
			if ((*chr_it).second[i]!=-1 && sampled_GC[(*chr_it).first][i]!=-1){
				int index = iround(sampled_GC[(*chr_it).first][i]);// (int)(sampled_GC[(*chr_it).first][i]*100);
				control_windows_count[index]++;
				density[index].push_back((*chr_it).second[i]);
		}
	}

	for(int i=0;i<101;i++)
		if(control_windows_count[i]>0){
			control_windows_count[i]=0;
			int ten_percent = density[i].size()/10;
			sort(density[i].begin(),density[i].end());
			float sum=0;
			for(unsigned int j=ten_percent;j<density[i].size()-ten_percent;j++){
				sum+=density[i][j];
				control_windows_count[i]++;
			}
			control_accum_density[i] = sum;
		}


//	startum_windows=0;
//	for(int i=0;i<21;i++){
//		GC_stratum+=control_accum_density[i];
//		startum_windows+=control_windows_count[i];
//	}


	for(int i=minGCPercentoConsider;i<=35;i+=2){
		GC_stratum=(control_accum_density[i]+control_accum_density[i-1])/(control_windows_count[i]+control_windows_count[i-1]);
		control_percent[i] = control_percent[i-1] = GC_stratum;
	}

	for(int i=0;i<minGCPercentoConsider;i++)
		control_percent[i] =control_percent[minGCPercentoConsider]; // = GC_stratum/startum_windows;


	for(int i=35;i<66;i++){
		control_percent[i]=control_accum_density[i]/control_windows_count[i];
	}


	for(int i = 66;i<maxGCPercentoConsider;i+=2){
		GC_stratum=(control_accum_density[i]+control_accum_density[i-1])/(control_windows_count[i]+control_windows_count[i-1]);
		control_percent[i] = control_percent[i-1] = GC_stratum;

	}


//	GC_stratum=0;
//	startum_windows=0;
//	for(int i=75;i<101;i++){
//		GC_stratum+=control_accum_density[i];
//		startum_windows+=control_windows_count[i];
//	}

	for(int i=maxGCPercentoConsider;i<101;i++)
		control_percent[i] = control_percent[maxGCPercentoConsider-1];// GC_stratum/startum_windows;
 }

	for (int i=0;i<101;i++){
		if (chip_percent[i]>max_chip)
			max_chip = chip_percent[i];

		 if (hasControl_ && control_percent[i]>max_control)
			max_control = control_percent[i];
	}

	int i ,chip_lower,chip_upper,control_lower,control_upper;
	i =0;
	while(chip_percent[i]<0.1*max_chip)
		i++;
	chip_lower =i;

	i=100;
	while(chip_percent[i]<0.1*max_chip)
		i--;
	chip_upper=i;

	if (hasControl_) {
        i=0;
        while(control_percent[i]<0.1*max_control)
            i++;
        control_lower =i;

        i=100;
        while(control_percent[i]<0.1*max_control)
            i--;
        control_upper=i;
	}

	int final_lower = chip_lower;
	int final_upper = chip_upper;

	if (hasControl_) {
        final_lower = max(control_lower,chip_lower);
        final_upper = min(control_upper,chip_upper);
	}

	int chip_total_windows=0,control_total_windows=0;
	for (int i=final_lower;i<=final_upper;i++){
		chip_total_windows+=chip_windows_count[i];
		if (hasControl_) control_total_windows+=control_windows_count[i];
	}


	// calculate the expectations
	for (int i=final_lower;i<=final_upper;i++){
		chip_lambda+=chip_percent[i]*chip_windows_count[i]/chip_total_windows;
		if (hasControl_) control_lambda+=control_percent[i]*control_windows_count[i]/control_total_windows;
	}

	float minScalingFactor = 0.25;
	float maxScalingFactor = 1/0.25;

	if (chip_lambda!=0 && !hasControl_ || chip_lambda!=0 && control_lambda!=0 && hasControl_){
//correcting here
		for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
            for(int i=0;i<sampled_sizes[(*chr_it).first];i++)
                if ((*chr_it).second[i] !=-1 && sampled_GC[(*chr_it).first][i]!=-1){
                        int index = (int)(sampled_GC[(*chr_it).first][i]*100);
                        //if (index>=final_lower && index<=final_upper) {  // will not check this condition starting from v1.43
							if (chip_percent[index]!=0) {
                                float scalingFactor=chip_lambda/chip_percent[index];
                                if (scalingFactor<minScalingFactor)
                                    scalingFactor=minScalingFactor;
                                if (scalingFactor>maxScalingFactor)
                                    scalingFactor=maxScalingFactor;
								(*chr_it).second[i] = (*chr_it).second[i]*scalingFactor;
							}
						//}
						//else
						//	(*chr_it).second[i]=-1;
				}
		}

        if (hasControl_)
            for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.end();++chr_it){
				for(int i=0;i< sampled_sizes[(*chr_it).first];i++)
					if ((*chr_it).second[i] !=-1 && sampled_GC[(*chr_it).first][i]!=-1){
						int index = (int)(sampled_GC[(*chr_it).first][i]*100);
//						if (index>=final_lower && index<=final_upper){ // will not check this condition starting from v1.43
							if (control_percent[index]!=0) {
                                float scalingFactor=control_lambda/control_percent[index];
                                if (scalingFactor<minScalingFactor)
                                    scalingFactor=minScalingFactor;
                                if (scalingFactor>maxScalingFactor)
                                    scalingFactor=maxScalingFactor;
								(*chr_it).second[i] = (*chr_it).second[i]*scalingFactor;
							}
//						}
						//else
						//	(*chr_it).second[i]=-1;
                    }
            }

        if (hasControl_ && calculateEmpiricalPvalue_) { //normalize also the controls to be used for the calculation of empirical p-values
            for(chr_it=sampled_control_density1.begin();chr_it!=sampled_control_density1.end();++chr_it){
				for(int i=0;i< sampled_sizes[(*chr_it).first];i++)
					if ((*chr_it).second[i] !=-1 && sampled_GC[(*chr_it).first][i]!=-1){
						int index = (int)(sampled_GC[(*chr_it).first][i]*100);
//						if (index>=final_lower && index<=final_upper){// will not check this condition starting from v1.43
							if (control_percent[index]!=0){
                                float scalingFactor=control_lambda/control_percent[index];
                                if (scalingFactor<minScalingFactor)
                                    scalingFactor=minScalingFactor;
                                if (scalingFactor>maxScalingFactor)
                                    scalingFactor=maxScalingFactor;
								(*chr_it).second[i] = (*chr_it).second[i]*scalingFactor;
							}
//						}
                    }
            }
            for(chr_it=sampled_control_density2.begin();chr_it!=sampled_control_density2.end();++chr_it){
				for(int i=0;i< sampled_sizes[(*chr_it).first];i++)
					if ((*chr_it).second[i] !=-1 && sampled_GC[(*chr_it).first][i]!=-1){
						int index = (int)(sampled_GC[(*chr_it).first][i]*100);
//						if (index>=final_lower && index<=final_upper){// will not check this condition starting from v1.43
							if (control_percent[index]!=0){
                                float scalingFactor=control_lambda/control_percent[index];
                                if (scalingFactor<minScalingFactor)
                                    scalingFactor=minScalingFactor;
                                if (scalingFactor>maxScalingFactor)
                                    scalingFactor=maxScalingFactor;
								(*chr_it).second[i] = (*chr_it).second[i]*scalingFactor;
							}
//						}
                    }
            }
        }



		// nosie ratio
		cout<<"chip lambda is: "<<chip_lambda<<endl;
		if (hasControl_) cout<<"control lambda is "<<control_lambda<<endl;

        if (hasControl_) {
            float noise_ratio = chip_lambda/control_lambda;
            if (noise_ratio>1)
                noise_ratio = 1;

            cout<<"noise ratio is: "<<noise_ratio<<endl;
            for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.end();++chr_it){
                for(int i=0;i< sampled_sizes[(*chr_it).first];i++)
                    if ((*chr_it).second[i]!=-1)
                        (*chr_it).second[i]=(*chr_it).second[i]*noise_ratio;

            }
        }
    } else {
		cerr<<"Warning: Sequencing Depth is too low can not perform GC content normalization and noise ratio adjustment. Skipping...."<<endl;
	}

}

void Profiler::calculate_GC_sampled_bins(string chr,string& chr_seq){
	sampled_GC[chr] = new float[sampled_sizes[chr]];

	for(int i=0;i<sampled_sizes[chr];i++){
		int position = i*bin_length; // up to version 1.43 (and I don't know why it was like this before): int position = i*bin_length+bin_length;
		int start = position-med < 0 ? 0:position-med;
		int length = position+med <= sizes[chr]? 2*med:sizes[chr]-position+1;
		float GC = calculate_GC(chr_seq,start,length);

		sampled_GC[chr][i] = GC;
	}

}



inline float Profiler::calculate_GC(string& chr_seq,int start,int length){
	//string window = chr_seq.substr(start,length);
	float G=0,C=0;

	if (length<=0){
		return -1;
	}
	for(int j =start;j<start+length;j++)
		if (chr_seq[j]!='N'){
			if (chr_seq[j] =='G' || chr_seq[j] == 'g')
				G++;
			if (chr_seq[j] == 'C' || chr_seq[j] == 'c')
				C++;
		}
		else{
			return -1;
		}

		return (G+C)/length;
}


void Profiler::derive_transition_probabilities_Input(){

	map<string,float *>::iterator chr_it;
	float transition_prob[2][2];
	transition_prob[0][0]=0;
	transition_prob[0][1]=0;
	transition_prob[1][0]=0;
	transition_prob[1][1]=0;
	for(chr_it=sampled_control_density1.begin();chr_it!=sampled_control_density1.end();++chr_it){
			int sampled_size = sampled_sizes[(*chr_it).first];
			for (int i=0;i<sampled_size-1;i++){
				if(peaks_Input[(*chr_it).first][i]==0 && peaks_Input[(*chr_it).first][i+1] == 0)
					transition_prob[0][0]+=1;
				else if(peaks_Input[(*chr_it).first][i]==0 && peaks_Input[(*chr_it).first][i+1] == 1)
					transition_prob[0][1]+=1;
				else if(peaks_Input[(*chr_it).first][i]==1 && peaks_Input[(*chr_it).first][i+1] == 0)
					transition_prob[1][0]+=1;
				else if(peaks_Input[(*chr_it).first][i]==1 && peaks_Input[(*chr_it).first][i+1] == 1)
					transition_prob[1][1]+=1;
			}

	}

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++){
			vector<float> row;
			row.push_back(i);
			row.push_back(j);
			row.push_back(transition_prob[i][j]);
			transition_Input.push_back(row);
		}
}

void Profiler::derive_transition_probabilities(){

	map<string,float *>::iterator chr_it;
	float transition_prob[2][2];
	transition_prob[0][0]=0;
	transition_prob[0][1]=0;
	transition_prob[1][0]=0;
	transition_prob[1][1]=0;
	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
			int sampled_size = sampled_sizes[(*chr_it).first];
			for (int i=0;i<sampled_size-1;i++){
				if(peaks[(*chr_it).first][i]==0 && peaks[(*chr_it).first][i+1] == 0)
					transition_prob[0][0]+=1;
				else if(peaks[(*chr_it).first][i]==0 && peaks[(*chr_it).first][i+1] == 1)
					transition_prob[0][1]+=1;
				else if(peaks[(*chr_it).first][i]==1 && peaks[(*chr_it).first][i+1] == 0)
					transition_prob[1][0]+=1;
				else if(peaks[(*chr_it).first][i]==1 && peaks[(*chr_it).first][i+1] == 1)
					transition_prob[1][1]+=1;
			}

	}

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++){
			vector<float> row;
			row.push_back(i);
			row.push_back(j);
			row.push_back(transition_prob[i][j]);
			transition.push_back(row);
		}
}





void Profiler::generate_observations(){
	map<string,float *>::iterator chr_it;
	observation_seq temp;


	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		temp.chr = (*chr_it).first;
		bool new_line=true;
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++){


			if (((*chr_it).second[i])!=-1 &&!hasControl_ || hasControl_ && ((*chr_it).second[i])!=-1 && sampled_control_density[(*chr_it).first][i]!=-1) {
				new_line = false;
				int obs;
				if (hasControl_)
                    obs = round(((*chr_it).second[i]-sampled_control_density[(*chr_it).first][i])); // difference
                else
                    obs = round(((*chr_it).second[i]));
				if (obs<0)
					temp.values.push_back(0);
				else
					temp.values.push_back(obs);
			} else {
				if (!new_line){
					temp.start = i-temp.values.size();
					new_line=true;
					if (temp.values.size()>0){
						observations.push_back(temp);
						temp.values.clear();
					}
				}
			}
        }
		temp.start = sampled_sizes[(*chr_it).first]-temp.values.size();
		new_line=true;
		if (temp.values.size()>0){
			observations.push_back(temp);
			temp.values.clear();

		}
	}
}

void Profiler::generate_observations_Input(){
	map<string,float *>::iterator chr_it;
	observation_seq temp;


	for(chr_it=sampled_control_density1.begin();chr_it!=sampled_control_density1.end();++chr_it){
		temp.chr = (*chr_it).first;
		bool new_line=true;
		for(int i=0;i<sampled_sizes[(*chr_it).first];i++){


			if (((*chr_it).second[i])!=-1 && sampled_control_density2[(*chr_it).first][i]!=-1) {
				new_line = false;
				int obs;
                obs = round(((*chr_it).second[i]-sampled_control_density2[(*chr_it).first][i])); // difference

				if (obs<0)
					temp.values.push_back(0);
				else
					temp.values.push_back(obs);
			} else {
				if (!new_line){
					temp.start = i-temp.values.size();
					new_line=true;
					if (temp.values.size()>0){
						observations_Input.push_back(temp);
						temp.values.clear();
					}
				}
			}
        }
		temp.start = sampled_sizes[(*chr_it).first]-temp.values.size();
		new_line=true;
		if (temp.values.size()>0){
			observations_Input.push_back(temp);
			temp.values.clear();

		}
	}
}

vector<observation_seq>  Profiler::get_obs(){
	return observations;
}
vector<observation_seq>  Profiler::get_obs_Input(){
	return observations_Input;
}



void Profiler::call_peaks(int mergeDist){

	map<string,float*>::iterator chr_it;

	int min_dist = mergeDist/bin_length+1;
    double lambda = 1;
	if (!hasControl_) {
        //	evaluate lambda: as round(mean density) of 70 or 90% of the highest density bins (ChIP and ATAC respectively)
        int maxNumberOfBinsToUse = 5000000;
        vector <float> densities;
        for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
            int sampled_size = sampled_sizes[(*chr_it).first];
            for(int i =0;i<sampled_size && densities.size()< maxNumberOfBinsToUse;i++)
                if(((*chr_it).second[i])>0){
                    densities.push_back((*chr_it).second[i]);
                }
        }
        float quantileToRemove = 1-ATACSEQPEAKPROP;
        if (isChIPseq_) quantileToRemove=1-CHIPSEQPEAKPROP;  // commented as it was removing too few "noise" peaks
        sort(densities.begin(), densities.end());
        double sum_of_elems = 0; int elCount=0;
        int maxElementToTake=floor(densities.size()*quantileToRemove);
        vector<float>::iterator nth = densities.begin() + maxElementToTake;
        for(std::vector<float>::iterator it = densities.begin(); it !=nth; ++it) {
            sum_of_elems += *it;elCount++;
        }
        if (elCount==0) {cerr << "Error: in the evaluation of lambda did not find any density values; please contact the HMCan developers\n"; exit(-1);}
        lambda =  sum_of_elems/elCount;
        densities.clear();
        cout << "Evaluated lambda for the Poisson test: "<<lambda<<endl;
        if (lambda==0) {
            cerr << "Error: cannot continue with lambda equal to zero\nAbort!\nPlease contact the HMCan developers"<<endl;
        }
	}


    for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
        int sampled_size = sampled_sizes[(*chr_it).first];
        for(int i =0;i<sampled_size;i++)
            peaks[(*chr_it).first][i]=0;
        for(int i =0;i<sampled_size;i++)
            if((*chr_it).second[i]!=-1){
                if (hasControl_) {
                    lambda = sampled_control_density[(*chr_it).first][i]>0 ? sampled_control_density[(*chr_it).first][i]:1;
                }
                float pvalue=poissoncdistribution(round((*chr_it).second[i])-1, lambda);
                if (pvalue<pvalue_threshold)
                    peaks[(*chr_it).first][i]=1;
            }
            else
                peaks[(*chr_it).first][i]=-1;
    }

	//remove single noise points
	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		int sampled_size = sampled_sizes[(*chr_it).first];
            for(int i =2;i<sampled_size-2;i++)
                if(peaks[(*chr_it).first][i]==1 && peaks[(*chr_it).first][i-1]==0  &&
							peaks[(*chr_it).first][i+1]==0)// && peaks[(*chr_it).first][i+2]==0)&& peaks[(*chr_it).first][i-2]==0
                    peaks[(*chr_it).first][i]=0;
	}

	for(chr_it=peaks.begin();chr_it!=peaks.end();++chr_it){
		int i=0;
		while(i<sampled_sizes[(*chr_it).first]-1){
			if((*chr_it).second[i]==1 && (*chr_it).second[i+1]==0){
				int count =0;
				int point =i+1;
				while((*chr_it).second[point]==0 && point<sampled_sizes[(*chr_it).first]){
					count++;
					point++;
				}
				if (count<=min_dist){
					for(int j=i+1;j<point;j++)
						peaks[(*chr_it).first][j]=1;
				}
				i=point;
			}
			else
				i++;
		}
	}
}


void Profiler::call_peaks_Input(int mergeDist){

	map<string,float*>::iterator chr_it;

	int min_dist = mergeDist/bin_length+1;
    double lambda = 1;


    for(chr_it=sampled_control_density1.begin();chr_it!=sampled_control_density1.end();++chr_it){
        int sampled_size = sampled_sizes[(*chr_it).first];
        for(int i =0;i<sampled_size;i++)
            peaks_Input[(*chr_it).first][i]=0;
        for(int i =0;i<sampled_size;i++)
            if((*chr_it).second[i]!=-1){
                lambda = sampled_control_density2[(*chr_it).first][i]>0 ? sampled_control_density2[(*chr_it).first][i]:1;

                float pvalue=poissoncdistribution(round((*chr_it).second[i])-1, lambda);
                if (pvalue<pvalue_threshold)
                    peaks_Input[(*chr_it).first][i]=1;
            }
            else
                peaks_Input[(*chr_it).first][i]=-1;
    }

	//remove single noise points
	for(chr_it=sampled_control_density1.begin();chr_it!=sampled_control_density1.end();++chr_it){
		int sampled_size = sampled_sizes[(*chr_it).first];
            for(int i =2;i<sampled_size-2;i++)
                if(peaks_Input[(*chr_it).first][i]==1 && peaks_Input[(*chr_it).first][i-1]==0  &&
							peaks_Input[(*chr_it).first][i+1]==0)
                    peaks_Input[(*chr_it).first][i]=0;
	}

	for(chr_it=peaks_Input.begin();chr_it!=peaks_Input.end();++chr_it){
		int i=0;
		while(i<sampled_sizes[(*chr_it).first]-1){
			if((*chr_it).second[i]==1 && (*chr_it).second[i+1]==0){
				int count =0;
				int point =i+1;
				while((*chr_it).second[point]==0 && point<sampled_sizes[(*chr_it).first]){
					count++;
					point++;
				}
				if (count<=min_dist){
					for(int j=i+1;j<point;j++)
						peaks_Input[(*chr_it).first][j]=1;
				}
				i=point;
			}
			else
				i++;
		}
	}
}


void Profiler::print_wig(string name){
    map<string,float*>::iterator chr_it;
    ofstream wig_file;
    ofstream pvalue_wig;
    wig_file.open((name+".wig").c_str());

    if(!wig_file){
    	cerr<<"Error: can not open WIG file for writing"<<endl;
    	exit(1);
    }
    wig_file<<"track name="<<"\""<<name<<"\" type=wiggle_0 visibility=2"<<endl;
    for(vector<observation_seq>::iterator ii=observations.begin();ii!=observations.end();++ii){
    	int start = (*ii).start*bin_length+1;
    	string chr = (*ii).chr;
    	wig_file<<"fixedStep chrom="<<chr<<" start="<<start<<" step="<<bin_length<<endl;

      //  cout << "DEBUG: "<<"fixedStep chrom="<<chr<<" start="<<start<<" step="<<bin_length<<endl;

        unsigned int numberOfElementsInObservation=(*ii).values.size();

     //   cout << "DEBUG: numberOfElementsInObservation = "<<numberOfElementsInObservation<<endl;
    //    cout << "DEBUG: maximal start position = "<<sizes[chr]-bin_length<<endl;

    	for(unsigned int i=0;i<numberOfElementsInObservation;i++){
    		if (start < sizes[chr]-bin_length)
    			wig_file<<(*ii).values[i]<<endl; //XXX
    		else
    			break;
    		start+=bin_length;
           // cerr << "DEBUG: printing current position: " << chr <<" : start " <<start <<endl;
    	}
    	wig_file<<endl;
    }
}


inline void Profiler::count_emissions(){
	int range,peak_index,obs_index;
	string chr;


	//restrict  observations to be in some specific range

	for (unsigned int i=0;i<observations.size();i++){
			for (unsigned int j=0;j<observations[i].values.size();j++){
				if (observations[i].values[j]<-10)
					observations[i].values[j] = -10;
			}
	}

	//get max and min obs
	for (unsigned int i=0;i<observations.size();i++){
		for (unsigned int j=0;j<observations[i].values.size();j++){
			if (observations[i].values[j]>max_obs   )
				max_obs = observations[i].values[j];
			if (observations[i].values[j]<min_obs)
				min_obs = observations[i].values[j];
		}
	}

	range = max_obs-min_obs+1;
	emission.push_back(vector<float>(range,0));
	emission.push_back(vector<float>(range,0));

	for (unsigned int i=0;i<observations.size();i++){
		chr = observations[i].chr;
		peak_index = observations[i].start;
		for (unsigned int j=0;j<observations[i].values.size();j++){
			obs_index = observations[i].values[j]-min_obs;
			if (peaks[chr][peak_index] == 1){

				if ( observations[i].values[j]>2)
					emission[1][obs_index]++;
				else
					emission[0][obs_index]++;
			}
			else if (peaks[chr][peak_index]==0){
				emission[0][obs_index]++;
			}
			peak_index++;
		}
	}
}


inline void Profiler::count_emissions_Input(){
	int range,peak_index,obs_index;
	string chr;


	//restrict  observations to be in some specific range

	for (unsigned int i=0;i<observations_Input.size();i++){
			for (unsigned int j=0;j<observations_Input[i].values.size();j++){
				if (observations_Input[i].values[j]<-10)
					observations_Input[i].values[j] = -10;
			}
	}

	//get max and min obs
	for (unsigned int i=0;i<observations_Input.size();i++){
		for (unsigned int j=0;j<observations_Input[i].values.size();j++){
			if (observations_Input[i].values[j]>max_obs_Input   )
				max_obs_Input = observations_Input[i].values[j];
			if (observations_Input[i].values[j]<min_obs_Input)
				min_obs_Input = observations_Input[i].values[j];
		}
	}

	range = max_obs_Input-min_obs_Input+1;
	emission_Input.push_back(vector<float>(range,0));
	emission_Input.push_back(vector<float>(range,0));

	for (unsigned int i=0;i<observations_Input.size();i++){
		chr = observations_Input[i].chr;
		peak_index = observations_Input[i].start;
		for (unsigned int j=0;j<observations_Input[i].values.size();j++){
			obs_index = observations_Input[i].values[j]-min_obs_Input;
			if (peaks_Input[chr][peak_index] == 1){

				if ( observations_Input[i].values[j]>2)
					emission_Input[1][obs_index]++;
				else
					emission_Input[0][obs_index]++;
			}
			else if (peaks_Input[chr][peak_index]==0){
				emission_Input[0][obs_index]++;
			}
			peak_index++;
		}
	}
}

int Profiler::get_minObs(){
	return min_obs;
}

int Profiler::get_minObs_Input(){
	return min_obs_Input;
}

vector<vector<float> > Profiler:: get_transition(){
	return transition;
}
vector<vector<float> > Profiler:: get_emission(){
	return emission;
}

vector<vector<float> > Profiler:: get_transition_Input(){
	return transition_Input;
}
vector<vector<float> > Profiler:: get_emission_Input(){
	return emission_Input;
}


void Profiler::clean_reads(std::map<std::string,std::vector<DNA_fragment > >& data,
			std::map<std::string,std::vector<DNA_fragment > >& control){

	map<string,vector<DNA_fragment > >::iterator chr_it;
	set<string> ToDelete;

	for (chr_it=data.begin();chr_it!=data.end();++chr_it){
		if(GC_profile.find((*chr_it).first) == GC_profile.end()){
			ToDelete.insert((*chr_it).first);
			cerr<<"Warning: "<<(*chr_it).first<<
					" can not be found on GC index, it will be ignored in further analysis"<<endl;
		}
	}

	for (chr_it=control.begin();chr_it!=control.end();++chr_it){
		if(GC_profile.find((*chr_it).first) == GC_profile.end()){
			ToDelete.insert((*chr_it).first);
			cerr<<"Warning: "<<(*chr_it).first<<
					" can not be found on GC index, it will be ignored in further analysis"<<endl;
		}
	}

	for (set<string>::iterator it=ToDelete.begin();it!=ToDelete.end();++it){
		data.erase((*it));
		control.erase((*it));
	}

	ToDelete.clear();


	vector<string> v1,v2,v3;

	for(chr_it=control.begin();chr_it!=control.end();++chr_it)
		v1.push_back((*chr_it).first);

	for (chr_it=data.begin();chr_it!=data.end();++chr_it)
		v2.push_back((*chr_it).first);

	sort(v1.begin(),v1.end());
	sort(v2.begin(),v2.end());


	v3 = vector<string>(v1.size()+v2.size(),"");
	vector<string>::iterator b = set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),v3.begin());

	v3.resize(b-v3.begin());


	for (chr_it=control.begin();chr_it!=control.end();++chr_it){
		vector<string>::iterator a = find(v3.begin(),v3.end(),(*chr_it).first);
		if (a==v3.end()){
			ToDelete.insert((*chr_it).first);
			cerr<<"warning: will delete "<<(*chr_it).first<<" since it is not common in all data"<<endl;
		}
	}

    if (control.begin()!=control.end()) {
        for (chr_it=data.begin();chr_it!=data.end();++chr_it){
            vector<string>::iterator a = find(v3.begin(),v3.end(),(*chr_it).first);
            if (a==v3.end()){
                ToDelete.insert((*chr_it).first);
                cerr<<"warning: will delete "<<(*chr_it).first<<" since it is not common in all data"<<endl;
            }
        }
	}


	for (set<string>::iterator it=ToDelete.begin();it!=ToDelete.end();++it){
			data.erase((*it));
			control.erase((*it));
		}



	return;
}



void Profiler::remove_blacklist(vector<BedEntry>& blacklist, std::map<std::string,float*>& profile){
	vector<BedEntry>::iterator chr_it;
	for (chr_it = blacklist.begin();chr_it!=blacklist.end();++chr_it){
		string chr = (*chr_it).chr;
		if(profile.find(chr) == profile.end())
			continue;
		int start_point = (*chr_it).start;
		while((start_point+1)%bin_length!=0)
			start_point++;

		for (int i= start_point;i<(*chr_it).end;i+=bin_length){
			int index = (i+1)/bin_length-1;
			if(index >= sampled_sizes[chr]){
				cerr<<"Blacklisted Value out of chromosome bounds. Exiting..."<<endl;
				exit(1);
			}
			profile[chr][index] = -1;

		}

	}
}

void Profiler::print_bedgraph(string name){

	string ChIP_bedgraph = name+"_ChIP.bedgraph";
	string Control_bedgraph = name+"_Input.bedgraph";
	ofstream ChIP,Control;
	map<string,float*>::iterator chr_it;


	ChIP.open(ChIP_bedgraph.c_str());
	Control.open(Control_bedgraph.c_str());


	if (!ChIP){
		cerr<<"Error: Can not open file: "<<ChIP_bedgraph<<". Exiting..."<<endl;
		exit(1);
	}

	if (!Control){
			cerr<<"Error: Can not open file: "<<ChIP_bedgraph<<". Exiting..."<<endl;
			exit(1);
	}

	ChIP<<"track type=bedGraph name=\""<<name+"_ChIP density"<<"\"description=\"Generated By HMCan\""<<
			"visibility=full color=200,100,0 altColor=0,100,200 priority=20"<<endl;

	Control<<"track type=bedGraph name=\""<<name+"_Input density"<<"\"description=\"Generated By HMCan\""<<
				"visibility=full color=200,100,0 altColor=0,100,200 priority=20"<<endl;

	for(chr_it=sampled_target_density.begin();chr_it!=sampled_target_density.end();++chr_it){
		string chr = (*chr_it).first;
		int sampled_size = sampled_sizes[chr];
		for(int i =0;i<sampled_size;i++)
			if (sampled_target_density[chr][i]!=-1)
				ChIP<<chr<<"\t"<<i*bin_length<<"\t"<<i*bin_length+bin_length-1<<"\t"<<sampled_target_density[chr][i]<<endl;
	}

	for(chr_it=sampled_control_density.begin();chr_it!=sampled_control_density.end();++chr_it){
			string chr = (*chr_it).first;
			int sampled_size = sampled_sizes[chr];
			for(int i =0;i<sampled_size;i++)
				if (sampled_control_density[chr][i]!=-1)
					Control<<chr<<"\t"<<i*bin_length<<"\t"<<i*bin_length+bin_length-1<<"\t"<<sampled_control_density[chr][i]<<endl;
		}

	ChIP.close();
	Control.close();

}

void Profiler::print_CNV_profile(string name){

    string CNV_profile_name = name+"_CNV_profile.txt";
    ofstream CNV_file;
    map<string,vector<float> >::iterator chr_it;

    CNV_file.open(CNV_profile_name.c_str());

    CNV_file<<"Chromosome\tStart\tRatio\tMedianRatio\tCopyNumber"<<endl;

    for (chr_it=medians.begin();chr_it!=medians.end();++chr_it){
        string chr_number = (*chr_it).first.substr(3); //skip chr prefix
        for (unsigned int i=0;i<(*chr_it).second.size();i++){
            int loci = i*large_bin_size+1;
            if ((*chr_it).second[i] == -1){// undefined bin
                CNV_file <<chr_number<<"\t"<<loci<<"\t-1\t-1\t-1"<<endl;
            }
            else{//define bin
                CNV_file<<chr_number<<"\t"<<loci<<"\t"<<ratio_profile[(*chr_it).first][i]
                        <<"\t"<<medians[(*chr_it).first][i]<<"\t-1"<<endl;
            }
        }


    }

    CNV_file.close();
}
