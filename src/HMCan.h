#ifndef HMCAN_H_
#define HMCAN_H_


#include <vector>
#include <string>
#include<cmath>
#include<cstdlib>
#include<fstream>

#include "utils.h"
#include "types.h"
class HMCan {
public:
	HMCan(std::vector<observation_seq>&);
	HMCan(std::vector<observation_seq>&,std::vector<std::vector<float> >&,std::vector<std::vector<float> >&,
			float ,float , int, int,int, float);
	virtual ~HMCan();
	void run(int,int);
	void addEmpericalPvalue (HMCan &, int fragment_length);
	void print(std::string name,int);
	void print_posterior(std::string name);

    std::vector<BedEntry> getFinalPeaks();
    std::vector<BedEntry> getFinalRegions();

private:

	int step,min_obs,bin_size;//, sample_size;
	float threshold1,threshold2,signal,background,posterior_threshold;
	std::vector<BedEntry> peaks,regions;
	std::vector<std::vector<float> > emissions,transitions;
	std::vector<int> states, samples_order;
	std::vector<observation_seq>& obs;

	void estimate_hmm_parameters(std::vector<std::vector<int> >&, std::vector<int>&,int );
	void merge_peaks(int);
	void print_probabilities(std::ofstream&,int);
	void forward_backward(std::vector<int>&, std::vector<double>&);
	int search_emission(int,int);
	int search_transition(int,int);
	void estimate_emissions();
	void posterior_decoding(std::vector<std::vector<int> >&, std::vector<int>&);
	void posterior_decoding_w_sampling(std::vector<std::vector<int> >&, std::vector<int>&);
	void estimate_transitions();
	void get_peaks(observation_seq&,std::vector<std::vector<int> >&,std::vector<int>&,int);
	void score_peak(BedEntry&,std::vector<int>&,std::vector<double>&);
	void distribution_info();
	float score_regions(int, int, int);
	void merge2();
	int myrandom (int );

};

#endif /* HMCAN_H_ */
