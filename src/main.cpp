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
#include <iostream>
#include <vector>
#include<iterator>
#include<algorithm>
#include <string>
#include <iomanip>
using namespace std;

#include "Profiler.h"
#include"Reader.h"
#include "HMCan.h"
#include "Parser.h"
const double HMCAN_VERSION = 1.44;

static void print_version()
{
    std::ostringstream ostr;
     ostr << "HMCan v" << std::fixed << std::setprecision(2) << HMCAN_VERSION << " : detection of histone modifications in normal and cancer ChIP-seq datasets\n";
     std::cout << ostr.str();
}


static void check_files(Parser p){
	ifstream test;

	//check GC index file
	test.open(p.GC_index.c_str());
	if (!test){
		cerr<<"Error: Can not find GC index file at "<<p.GC_index<<endl;
		exit(1);
	}
}


int main(int argc, char * argv[]){

	print_version();
	if (argc<5){
		cout<<"Not enough arguments"<<endl<<"Usage: ./HMCan <TargetFile> <ControlFile> <configuration file> <Name>"<<endl;
        cout<<"\nIf you do not have a control file: type \"-\", e.g.,\n    ./HMCan <TargetFile> - <configuration file> <Name>"<<endl;
		cout<<"For more details please see README.md"<<endl;
		return 1;
	}
	Parser p = Parser();
	p.parse(argv[3]);
	string name = argv[4];

	check_files(p);

    vector<vector<float> > emission,transition,emission_Input,transition_Input;
	vector<observation_seq> obs, obs_Input;

	bool hasControl=true;
    if (strcmp(argv[2],"-")==0 || !is_file_exist(argv[2])) {
        hasControl=false;
        cout << "Warning: HMCan is not going to use a control dataset for this analysis\n";
    }

    Profiler test_gc = 	Profiler(p.GC_index.c_str(),p.min_length,p.median_length,p.max_length ,p.bin_size,p.large_bin_size,p.merge_dist,p.pvalue_threshold,p.pairedEnds,p.isCHIPseq,hasControl, p.GC_merge_dist, p.calculateEmpiricalPvalue);

    if (string(argv[1]).substr( string(argv[1]).length() - 3 ) == "wig" && p.t!=WIG) {
        cerr << "Warning: the type of input file was set to .wig\n";
        p.t=WIG;
    }

    if (string(argv[1]).substr( string(argv[1]).length() - 3 ) == "bed" && p.t!=BED) {
        cerr << "Warning: the type of input file was set to .bed\n";
        p.t=BED;
    }

    if (string(argv[1]).substr( string(argv[1]).length() - 3 ) == "bam" && p.t!=BAM) {
        cerr << "Warning: the type of input file was set to .bam\n";
        p.t=BAM;
    }

    if (string(argv[1]).substr( string(argv[1]).length() - 3 ) == "sam" && p.t!=SAM) {
        cerr << "Warning: the type of input file was set to .sam\n";
        p.t=SAM;
    }

	if (p.t != WIG) { //read BAM, SAM or BED files:
        Reader data(p.t,10,p.removeDuplicates);    //minimal read quality of 10
        Reader control(p.t,10,p.removeDuplicates); //minimal read quality of 10
        cout<<"Reading chip...."<<endl;
        int readsRead = 0;
        if (!is_file_exist(argv[1])) {
            cerr << "Error: the input file "<<argv[1]<<" does not exist\nCheck your command line and rerun\n";
            exit (-1);
        }
        map<string,vector<DNA_fragment> > target_data = data.Read(argv[1],readsRead);
        cout<<"Kept "<<readsRead<<" reads"<<endl;
        map<string,vector<DNA_fragment> > control_data ;

        if (hasControl) {
            cout<<"Reading control...."<<endl;
            control_data = control.Read(argv[2],readsRead);
            cout<<"Kept "<<readsRead<<" reads"<<endl;
        }
        vector<BedEntry> black_list = data.Read_blacklist(p.blacklist.c_str());
        test_gc.build_profile(target_data,control_data,p.path,black_list,p.ifKeepChIPdensityIntact,p.CNAnormalization);
        target_data.clear();
        control_data.clear();
    } else { //read density from the provided wig file
    	test_gc.build_profile_wig(p.genomeLengthFile.c_str(),argv[1]);
    }

	int step = (int) ceil(float(test_gc.getMaxLength())/p.bin_size);

	if (p.calculateEmpiricalPvalue && p.callpeaks && p.t != WIG && hasControl) {
        emission_Input = test_gc.get_emission_Input();
        transition_Input = test_gc.get_transition_Input();
        obs_Input = test_gc.get_obs_Input();
        HMCan i = HMCan(obs_Input,transition_Input,emission_Input,p.iter_threshold,p.final_threshold,step,test_gc.get_minObs_Input(),p.bin_size,p.posterior_threshold);
        emission = test_gc.get_emission();
        transition = test_gc.get_transition();
        obs = test_gc.get_obs();
        HMCan a = HMCan(obs,transition,emission,p.iter_threshold,p.final_threshold,step,test_gc.get_minObs(),p.bin_size,p.posterior_threshold);
        a.run(p.merge_dist,p.max_iter);
        cout << "Running HMM for the Input mock data\n";
        i.run(p.merge_dist,p.max_iter);

        a.addEmpericalPvalue(i, p.median_length);

        a.print(name,p.median_length);
        if (p.posterior){
            a.print_posterior(name);
        }

	} else if (p.callpeaks){

        emission = test_gc.get_emission();
        transition = test_gc.get_transition();
        obs = test_gc.get_obs();
        HMCan a = HMCan(obs,transition,emission,p.iter_threshold,p.final_threshold,step,test_gc.get_minObs(),p.bin_size,p.posterior_threshold);

        a.run(p.merge_dist,p.max_iter);
        a.print(name,p.median_length);
        if (p.posterior){
            a.print_posterior(name);
        }

	}

	if (p.wig && p.t != WIG)
		test_gc.print_wig(name);

	if(p.bedgraph && p.t != WIG)
		test_gc.print_bedgraph(name);
    if (p.CNAnormalization)
        test_gc.print_CNV_profile(name);


	return 0;

}
