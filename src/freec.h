#ifndef FREEC_H_
#define FREEC_H_

#include<string>
#include<vector>
#include<map>

void calculateCopyNumberMedians(int ploidy,std::map <std::string,
								std::vector <int> > &breakpoints, std::map <std::string, std::vector <float> > &medians,
								std::map <std::string, std::vector <float> > &ratio_profile);



int recalculateRatioUsingCG (int ploidy, std::map <std::string, std::vector <float> > &read_count,std::map <std::string,
							std::vector <float> > &GC_profile,std::map <std::string, std::vector <float> > &notNprofile,
							std::map <std::string, std::vector <float> > &ratio_profile);


double calculateMedianAround (float interval, float around, std::map <std::string, std::vector <float> > &read_count,
							 std::map <std::string, std::vector <float> > &GC_profile );


void calculateBreakpoints(std::map <std::string, std::vector <float> > &ratio_profile,
		std::map <std::string, std::vector <int> > &breakpoints, double breakpointThreshold);

void removeOutliers (std::vector <float> &ratio_profile, float RangeMin, float RangeMax) ;



#endif /* FREEC_H_ */
