#ifndef PROFILER_H_
#define PROFILER_H_

#define CHIPSEQPEAKPROP 0.3
#define ATACSEQPEAKPROP 0.1


#include <map>
#include <string>
#include<vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <numeric>
#include <iostream>
#include <cstdlib>
#include "utils.h"

class Profiler {
public:
	Profiler(const char* ,int , int, int ,int ,int ,int,float,bool, bool, bool hasControl, int, bool);
	virtual ~Profiler();
	void build_profile(std::map<std::string,std::vector<DNA_fragment > >& , std::map<std::string,std::vector<DNA_fragment > >&,
			std::string&, std::vector<BedEntry>&, bool, bool);
    void build_profile_wig(const char *, const char *);
	std::vector<std::vector<float> > get_transition();
	std::vector<std::vector<float> > get_emission();
	std::vector<observation_seq>  get_obs();
    std::vector<std::vector<float> > get_transition_Input();
	std::vector<std::vector<float> > get_emission_Input();
	std::vector<observation_seq>  get_obs_Input();

	int get_minObs();
    int get_minObs_Input();

	void print_wig(std::string);
	void print_bedgraph(std::string);
    	int getMaxLength();
    	void calculateFragLengthDist(std::vector<DNA_fragment> &);
    	void print_CNV_profile(std::string);
private:

	int med,bin_length,left,right,large_bin_size,merge_dist,min_obs,max_obs, GC_mergeDist,min_obs_Input, max_obs_Input;
	float h,slope1,slope2,b1,b2,left_area,right_area,reads_ratio,pvalue_threshold;
	bool pairedEnds_,isChIPseq_,hasControl_,calculateEmpiricalPvalue_;
	std::vector<std::vector<float> > transition,emission, transition_Input,emission_Input;
	std::map<std::string,std::vector<float> > GC_profile,medians,notNprofile,ratio_profile;
	std::map<std::string,int > sizes,sampled_sizes;
	std::map<std::string,float*> sampled_target_density,sampled_control_density,sampled_GC,peaks,sampled_control_density1,sampled_control_density2, peaks_Input;
	std::vector<observation_seq> observations,observations_Input;



	void build_single_profile(std::vector<DNA_fragment >& , std::string, bool, bool);
	void read_gc_profile(std::ifstream& );
	void calculate_nonNs_ratio(std::string chr, std::string& chr_seq);
	void call_freec(std::map<std::string,std::vector<DNA_fragment > >&, bool isInput);
	void extend_reads(std::vector<DNA_fragment >&,int);
	int getMaxFragLength(std::vector<DNA_fragment >&);
	void calculate_density_coefs();
	void Normalize_GC();
	void calculate_GC_sampled_bins(std::string,std::string&);
	inline float calculate_GC(std::string&,int,int);
	void print_values();
	void print_vector(float *, int, const char *);
	void derive_transition_probabilities();
    void derive_transition_probabilities_Input();

	//void derive_emision_probabilities();
	inline void count_emissions();
	inline void count_emissions_Input();


	void generate_observations();
    void generate_observations_Input();

	void call_peaks(int mergeDist);
    void call_peaks_Input(int mergeDist);

	void peaks_ideal();
	void clean_reads(std::map<std::string,std::vector<DNA_fragment > >&,
			std::map<std::string,std::vector<DNA_fragment > >&);
	void remove_blacklist(std::vector<BedEntry>&,std::map<std::string,float*>& );
    void getFragmentLengthDistr (std::vector<int> &) ;


};

#endif /* PROFILER_H_ */
