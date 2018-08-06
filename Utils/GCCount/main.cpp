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


// main.cpp

#include "SVfinder.h"
#include "GenomeCopyNumber.h"

using namespace std ;
int verbose = false;

int main(int argc, char *argv[]) {
	//check arguments
    cout << "gccount: calculate GC-profile\n";
	if (argc < 3) {
		cerr << "\n\tPlease specify a config file:\n./gccount -conf config.txt\n";
		return 0;
	}
	if (argc == 3) {
		if (
			(strcmp(argv[1], "-conf") != 0)&&
			(strcmp(argv[1], "--conf") != 0)
			) {
            cerr << "\n\tPlease specify a config file\n\n";
			return 0;
		}
	}
	//check if config exists
    ifstream ifile(argv[2]);
    if (!ifile) {
      cerr << "\n\tCould not find your config file.. Please, check the existance of "<<argv[2] <<"\n\n";
      return -1;
    }

	//read config file
	ConfigFile cf(argv[2]);

    int window;

	try {
		window = int(cf.Value("general","window"));
        cout << "..windows size was set to "<< window<< "\n";

	} catch (const char * error) {
        cerr << "you need to specify window size\n";
        exit(0);
	}


	string myName = "";
	string controlName = "";
	string outputDir = ".";


	try {
		outputDir = std::string(cf.Value("general","outputDir"));
        if ( access( outputDir.c_str(), 0 ) == 0 )
        {
            struct stat status;
            stat( outputDir.c_str(), &status );

            if ( status.st_mode & S_IFDIR )
            {
                cout << "..Output directory:\t" << outputDir << "\n";
            }
            else
            {
                cerr << "Error: The path you entered for 'outputDir': "<< outputDir <<" is a file. It shoud be a directory" << endl;
                exit(-1);
            }
        }
        else
        {
            cerr << "Error: Path "<<outputDir<<" doesn't exist." << endl;
            exit(-1);
        }


	} catch (const char * error) {
		//Do nothing, outputDir will be "."
	}


    string dirWithFastaSeq = "";

    try {
            dirWithFastaSeq = std::string(cf.Value("general","chrFiles"));

            if ( access( dirWithFastaSeq.c_str(), 0 ) == 0 )
            {
                struct stat status;
                stat( dirWithFastaSeq.c_str(), &status );

                if ( status.st_mode & S_IFDIR )
                {
                    cout << "..Directory with files containing chromosome sequences:\t" << dirWithFastaSeq << "\n";
                }
                else
                {
                    cerr << "Error: The path you entered for 'dirWithFastaSeq': "<< dirWithFastaSeq <<" is a file. It shoud be a directory" << endl;
                    exit(-1);
                }
            }
            else
            {
                cerr << "Error: Path "<<dirWithFastaSeq<<" doesn't exist. Comment the line with 'chrFiles' if you use a precalculated GC-content profile or set the correct path" << endl;
                exit(-1);
            }
    } catch (const char * error) {
            cerr << "..you need to provide a directory with fasta files. Unable to proceed\n";
            exit(-1);
    }



    int step = NA;
    try {
        step =  (int)cf.Value("general","step");
		cout << "..Step:\t" << step<< "\n";
	} catch (const char * error) {
		//Do nothing
	}

    //READ SAMPLE DATA:

	GenomeCopyNumber sampleCopyNumber;

	if (step != NA) {
        sampleCopyNumber.setStep(step);
	}


    try {
            cout << "..File with chromosome lengths:\t" << std::string(cf.Value("general","chrLenFile")) << "\n";
    } catch (const char * error) {
        cerr << "You need to provide a file with chromosome lengths\n";
        exit (-1);
    }

    try {
        sampleCopyNumber.readChrInfo(cf.Value("general","chrLenFile"),window, step);
    }catch (const char * error) {
        cerr << "Something went wrong :(\n";
        exit (-1);
    }

	if (step == NA) {
        step= window;
	}

    string GCprofileFile = "";

	sampleCopyNumber.fillCGprofile(dirWithFastaSeq);
	GCprofileFile = outputDir+"/GC_profile.cnp";

    try { //read mappability file
        string mapFile = cf.Value("general","gemMappabilityFile");
        sampleCopyNumber.readGemMappabilityFile(mapFile);
            //rewrite GC-profile with mappability as the last (5th) colomn

        cout << "..Mappability track from "<<mapFile<<" has been added to "<< GCprofileFile <<"\n";
    } catch (const char * error) {
            //
    }
	sampleCopyNumber.printCGprofile(GCprofileFile);
	return 0;
}
