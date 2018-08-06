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

using namespace std;

#include "freec.h"
#include "segmentation.h"

//map <string, vector <float> > read_count;
//map <string, vector <float> > ratio_profile;
//map <string, vector <float> > GC_profile;
//map <string, vector <float> > medians;
//map <string, vector <float> > notNprofile;//percentage of not 'N' in a window
//map <string, vector <int> > breakpoints;



void calculateCopyNumberMedians(int ploidy,std::map <string, vector <int> > &breakpoints,
		std::map <string, vector <float> > &medians,std::map <string, vector <float> > &ratio_profile) {
    map<string,vector <float> >::iterator it;
	for ( it=ratio_profile.begin() ; it != ratio_profile.end(); it++ ) {
        int breakPointStart = 0;
        int breakPointEnd;
        float median;
        int length_ =  (*it).second.size();
        string chrom = (*it).first;
        breakpoints[chrom].push_back(length_-1);

        vector <float>medianProfile_ (length_,-1);

        for (int i  = 0; i < (int)breakpoints[chrom].size();i++) {
            breakPointEnd =  breakpoints[chrom][i];
            vector<float> data;
            int notNA = 0;
            for (int j = breakPointStart; j <= breakPointEnd; j++)
                if (ratio_profile[chrom][j] != -1) {
                    data.push_back(ratio_profile[chrom][j]);
                    notNA++;
                }
            int totalCount = breakPointEnd-breakPointStart+1;
            if(data.size()==0){
            	median = NA;

            } else
            	median = get_median(data); //including the last point



            for (int j = breakPointStart; j<= breakPointEnd; j++) {
                medianProfile_[j] = median;
            }

            breakPointStart = breakPointEnd+1;
            data.clear();
        }
        medians[chrom] = medianProfile_;
        medianProfile_.clear();
	}
}

int recalculateRatioUsingCG (int ploidy, map <string, vector <float> > &read_count,map <string, vector <float> > &GC_profile,
							  map <string, vector <float> > &notNprofile,map <string, vector <float> > &ratio_profile) {


    bool intercept = 1;
    float minExpectedGC = .35;
    float maxExpectedGC = .55;
    float minMappabilityPerWindow = 0.6;


    //try degree 3 and 4, SELECT THE ONE WITH LESS ITTERATIONS:
    int minDegreeToTest = 3;
    int maxDegreeToTest = 4;
    int selectedDegree=minDegreeToTest;
    int maximalNumberOfIterations = 200;
    int bestNumberOfIterations=maximalNumberOfIterations;
    int maximalNumberOfCopies = ploidy*2;
    float interval = float (0.01) ;


    vector <float> y; //y ~ ax^2+bx+c
    vector <float> x;
    //fill x and y:
    map<string,vector <float> >::iterator it;
    for ( it=read_count.begin() ; it != read_count.end(); it++ ) {
        if (! (((*it).first).find("X")!=string::npos || ((*it).first).find("Y")!=string::npos)){ //??
		for (int i = 0; i< (*it).second.size(); i++) {
               		float value = (*it).second.at(i);
                	float nonN = notNprofile[(*it).first].at(i);
                	float GC = GC_profile[(*it).first].at(i);
                	if ((value>0)&&(nonN>minMappabilityPerWindow)) {
                    		x.push_back(GC);
                    		y.push_back(value/nonN);
                	}
		}
        }
     }

    for (int degree=minDegreeToTest;degree<=maxDegreeToTest; degree++) {

        //first guess about parameters
        const int npoints = degree+2;

        double around [MAXDEGREE+2];
        for (int i = 0; i<npoints; i++) {
            around[i] = minExpectedGC + (maxExpectedGC-minExpectedGC)/(npoints-1)*i; //0.55-0.35
        }  // for degree = 3 one will get : { 0.35, 0.40, 0.45, 0.5, 0.55 }; for 0.35 and 0.55


        double yValues [MAXDEGREE+2];
        for (int i = 0; i <npoints; i++)
            yValues[i] = calculateMedianAround(interval, float(around[i]),read_count,GC_profile);
        int nvars = degree; //fit by cubic polynomial
        ap::real_2d_array xy;
        xy.setlength(npoints,nvars+1);
        for (int i = 0; i <npoints; i++) {
            xy(i,degree) = yValues[i];
            xy(i,degree-1) = around[i];
            for (int j = degree-2; j>=0; j--) {
                xy(i,j)=xy(i,j+1)*around[i];
            }
            /* this is equal to
            xy(i,0) = around[i]*around[i]*around[i];
            xy(i,1) = around[i]*around[i];
            xy(i,2) = around[i];
            xy(i,3) = yValues[i]; */
        }
        linearmodel lm;
        int info;
        lrreport ar;

        if  (intercept)
            lrbuild(xy,npoints,nvars,info,lm,ar);
        else
            lrbuildz(xy,npoints,nvars,info,lm,ar);
        if (info != 1) {
            cerr << "Error in the first linear regression (the first guess about parameters), code: " << info <<"\n";
        }
        ap::real_1d_array v;
        v.setlength(nvars+int(intercept));
        lrunpack(lm,v,nvars);
        double a[MAXDEGREE+1];
        for (int i = 0; i <degree; i++) {
                a[i] = v(i);
        }

        if  (intercept)
            a[degree] = v(degree);
        else
            a[degree] = 0;

        int realNumberOfIterations = maximalNumberOfIterations;
        float rmserror = runEM(x,y,a,degree,realNumberOfIterations,ploidy,maximalNumberOfCopies, intercept);
        if (rmserror == -1) {
            cerr << "Error in EM => unable to calculate normalized profile\n";
            return NA;
        }
        //cout << "root mean square error = " << rmserror << "\n";
        if (realNumberOfIterations<bestNumberOfIterations) {
                selectedDegree = degree;
                bestNumberOfIterations=realNumberOfIterations;
        }
        cout << "Y = ";
        for (int i=0; i<degree;i++) {
            cout << a[i]<<"*x^" <<  degree-i <<"+" ;
        }
        cout << a[degree] <<"\n";

    }


    int degree = selectedDegree;
    const int npoints = degree+2;

    double around [MAXDEGREE+2];
    for (int i = 0; i<npoints; i++) {
            around[i] = minExpectedGC + (maxExpectedGC-minExpectedGC)/(npoints-1)*i; //0.55-0.35
    }  // for degree = 3 one will get : { 0.35, 0.40, 0.45, 0.5, 0.55 }; for 0.35 and 0.55


    double yValues [MAXDEGREE+2];
    for (int i = 0; i <npoints; i++)
        yValues[i] = calculateMedianAround(interval, float(around[i]),read_count,GC_profile);
    int nvars = degree; //fit by cubic polynomial
    ap::real_2d_array xy;
    xy.setlength(npoints,nvars+1);
    for (int i = 0; i <npoints; i++) {
            xy(i,degree) = yValues[i];
            xy(i,degree-1) = around[i];
            for (int j = degree-2; j>=0; j--) {
                xy(i,j)=xy(i,j+1)*around[i];
            }
            /* this is equal to
            xy(i,0) = around[i]*around[i]*around[i];
            xy(i,1) = around[i]*around[i];
            xy(i,2) = around[i];
            xy(i,3) = yValues[i]; */
    }
    linearmodel lm;
    int info;
    lrreport ar;

    if  (intercept)
            lrbuild(xy,npoints,nvars,info,lm,ar);
    else
            lrbuildz(xy,npoints,nvars,info,lm,ar);
    if (info != 1) {
        cerr << "Error in the first linear regression (the first guess about parameters), code: " << info <<"\n";
    }
    ap::real_1d_array v;
    v.setlength(nvars+int(intercept));
    lrunpack(lm,v,nvars);
    double a[MAXDEGREE+1];
    for (int i = 0; i <degree; i++) {
                a[i] = v(i);
    }

    if  (intercept)
            a[degree] = v(degree);
    else
            a[degree] = 0;

    bestNumberOfIterations = maximalNumberOfIterations;
    float rmserror = runEM(x,y,a,degree,bestNumberOfIterations,ploidy,maximalNumberOfCopies, intercept);
    if (rmserror == -1) {
            cerr << "Error in EM => unable to calculate normalized profile\n";
            return NA;
    }
	//cout << "root mean square error = " << rmserror << "\n";
    cout << "Y = ";
    for (int i=0; i<degree;i++) {
       	cout << a[i]<<"*x^" <<  degree-i <<"+" ;
	}

    cout << a[degree] <<"\n";

    //update v1.16:

    // get boundaries on GC-content where the normalization will work well

	vector <float> res (maximalNumberOfCopies+1);
    vector<vector<float> > residue(101);

    for ( it=read_count.begin() ; it != read_count.end(); it++ ) {
        for (int i = 0; i< (*it).second.size(); i++) {
            float gc = GC_profile[(*it).first].at(i);
            float value = (*it).second.at(i);
            float nonN = notNprofile[(*it).first].at(i);
            if ((gc>0)&&(nonN>minMappabilityPerWindow) && (value >=0)) {
                for (int j = 0; j <= maximalNumberOfCopies; j++)
                    res[j] = fabs(polynomial(gc,a,float(j)/ploidy,degree)-value);
                float min_resed = res[get_min_index(res)];
                int intGC = round(gc*100);
                residue[intGC].push_back(min_resed);
            }
        }
    }
    vector <float> MeanDeviationsFromTheFit;
    float thresholdOfDeviationFromTheFit;
    for (int intGC = 1;intGC<100;intGC++) {
        if (residue[intGC].size()>0) {
            MeanDeviationsFromTheFit.push_back(get_mean(residue[intGC]));
        }
    }
    float upperQuartile = get_upper_quartile(MeanDeviationsFromTheFit);
    float IQR = get_iqr(MeanDeviationsFromTheFit);
    thresholdOfDeviationFromTheFit = 1.5*IQR+upperQuartile;
    cout << "IQR: "<<IQR<<"; upperQ: "<<upperQuartile<<"; threshold: "<<thresholdOfDeviationFromTheFit<<endl;
    float minGCtoConsider = 0;
    float maxGCtoConsider = 1;

    for (int intGC = 50;intGC<100;intGC++) {
        if (residue[intGC].size()>0) {
            if (get_mean(residue[intGC])>thresholdOfDeviationFromTheFit) {
                maxGCtoConsider = intGC*1.0/100;
                break;
            }
        }
    }
    for (int intGC = 45;intGC>0;intGC--) {
        if (residue[intGC].size()>0) {
            if (get_mean(residue[intGC])>thresholdOfDeviationFromTheFit) {
                minGCtoConsider = intGC*1.0/100;
                break;
            }
        }
    }
    cout << "HMCan will only calculate ratio for genomic windows with GC-content from "<< minGCtoConsider << " to "<<maxGCtoConsider<<endl;
    //end of update

    for ( it=read_count.begin() ; it != read_count.end(); it++ ) {
    	cout<<(*it).first<<'\t'<<GC_profile[(*it).first].size()<<'\t'<<notNprofile[(*it).first].size()<<endl;
        if ((int)ratio_profile[(*it).first].size()!=((*it).second).size())
            ratio_profile[(*it).first].resize(((*it).second).size());
        for (int i = 0; i< (*it).second.size(); i++) {
            float x = GC_profile[(*it).first].at(i);
            float value = (*it).second.at(i);
            float nonN = notNprofile[(*it).first].at(i);
            float rati;
            //if ((x>0)&&(nonN>minMappabilityPerWindow))
            if ((x>=minGCtoConsider)&&(x<=maxGCtoConsider)&&(nonN>minMappabilityPerWindow) && value>=0) //update v1.16
                rati =value/polynomial(x,a,1,degree)/nonN;
            else{

                rati = NA;
            }
            if ((rati != NA)&&(rati < 0))
                rati = 0; //this happens if  polynomial(x,a,b,c,d) < 0
            ratio_profile[(*it).first].at(i) =rati ;
        }
    }
    x.clear();
    y.clear();
    return bestNumberOfIterations;
}

void removeOutliers (std::vector <float> &ratio_profile, float RangeMin, float RangeMax) {
    if (ratio_profile.size()<10)
        return;
    for (int i = 1; i< (ratio_profile.size()-1); i++) {
        if (ratio_profile[i]>RangeMin*max(ratio_profile[i-1],ratio_profile[i+1]) && ratio_profile[i]<RangeMax*max(ratio_profile[i-1],ratio_profile[i+1])) {
            ratio_profile[i]=NA;
            i++; //skip the next element
        }
    }
}

double calculateMedianAround (float interval, float around, map <string, vector <float> > &read_count,
							 map <string, vector <float> > &GC_profile ) {

	float maxCG = around+interval;
	float minCG = around-interval;

	vector <float> myValuesAround;
    map<string,vector <float> >::iterator it;
	for ( it=read_count.begin() ; it != read_count.end(); it++ ) {
        if (! (((*it).first).find("X")!=string::npos || ((*it).first).find("Y")!=string::npos))
			for (int i = 0; i< (*it).second.size(); i++) {
                float GC = GC_profile[(*it).first].at(i);
				if ((  GC<=maxCG)&&(GC>=minCG)) {
                    float value = (*it).second.at(i);
					if (value>0) //non-zero values
						myValuesAround.push_back(value);
                }
			}
	}
    if (myValuesAround.size()==0) {
        cerr << "Error: zero reads in windows with the GC-content around "<<around<< " with interval "<< interval<<", will try again with "<< interval*4<<"\n";
        cerr << "Unable to proceed..\n";
        cerr << "Try to rerun the program with higher number of reads\n";
        exit(-1);

    }
	float median = get_median(myValuesAround);
	myValuesAround.clear();
	return median;
}

void calculateBreakpoints(std::map <string, vector <float> > &ratio_profile, std::map <string, vector <int> > &breakpoints, double breakPointThreshold) {
    int breakPointType = 2;
	cout << "..Calculating breakpoints, breakPointThreshold = " <<breakPointThreshold<<"\n";
	map<string,vector <float> >::iterator it;
	for ( it=ratio_profile.begin() ; it != ratio_profile.end(); it++ ) {
        cout << "..processing chromosome " <<(*it).first<<"\n";
        int result = calculateBreakpoints_general(breakPointThreshold,(*it).second.size(),(*it).second,breakpoints[(*it).first],0,breakPointType);
		if (result == 0) {
            cerr << "..failed to run segmentation on chr" << (*it).first << "\n";
		}
	}
}
