HMCan: main.o Reader.o Profiler.o  ap.o igammaf.o gammafunc.o normaldistr.o ablasf.o ablas.o ortfac.o sblas.o rotations.o reflections.o linreg.o hblas.o descriptivestatistics.o creflections.o blas.o bdsvd.o svd.o ialglib.o utils.o segmentation.o freec.o poissondistr.o HMCan.o Parser.o
	g++ -g -O3  -Wall  -m64  -lm -o HMCan main.o ap.o igammaf.o gammafunc.o normaldistr.o ablasf.o ablas.o ortfac.o sblas.o rotations.o reflections.o linreg.o hblas.o descriptivestatistics.o creflections.o blas.o bdsvd.o svd.o ialglib.o Profiler.o Reader.o utils.o segmentation.o freec.o poissondistr.o HMCan.o Parser.o
	rm *.o
main.o: main.cpp 
	g++ -O3 -g  -c -Wall -m64 main.cpp
ap.o: ap.cpp ap.h
	g++ -c -O3 -g  -Wall -m64  ap.cpp
igammaf.o: igammaf.cpp igammaf.h
	g++ -c -O3 -g  -Wall -m64  igammaf.cpp
gammafunc.o: gammafunc.cpp gammafunc.h
	g++ -c -O3 -g   -Wall -m64  gammafunc.cpp
normaldistr.o: normaldistr.cpp normaldistr.h
	g++ -c -O3 -g  -Wall -m64  normaldistr.cpp
ablasf.o: ablasf.cpp ablasf.h
	g++ -c -O3 -g  -Wall -m64  ablasf.cpp
ablas.o: ablas.cpp ablas.h
	g++  -c -O3 -g  -Wall -m64  ablas.cpp
ortfac.o: ortfac.cpp ortfac.h
	g++  -c -O3 -g  -Wall -m64  ortfac.cpp
sblas.o: sblas.cpp sblas.h
	g++  -c -O3 -g  -Wall -m64  sblas.cpp
rotations.o: rotations.cpp rotations.h
	g++  -c -O3 -g  -Wall -m64  rotations.cpp
reflections.o:  reflections.cpp reflections.h
	g++  -c -O3 -g  -Wall -m64  reflections.cpp
linreg.o:  linreg.cpp linreg.h
	g++  -c -O3 -g  -Wall -m64 linreg.cpp
hblas.o:  hblas.cpp hblas.h
	g++  -c -O3 -g  -Wall -m64  hblas.cpp
descriptivestatistics.o:  descriptivestatistics.cpp descriptivestatistics.h
	g++  -c -O3 -g  -Wall -m64 descriptivestatistics.cpp 
creflections.o:  creflections.cpp creflections.h
	g++  -c -O3 -g  -Wall -m64 creflections.cpp
blas.o:  blas.cpp blas.h
	g++  -c -O3 -g  -Wall -m64 blas.cpp
bdsvd.o:  bdsvd.cpp bdsvd.h
	g++  -c -O3 -g  -Wall -m64 bdsvd.cpp
svd.o:  svd.cpp svd.h
	g++  -c -O3 -g  -Wall -m64 svd.cpp
ialglib.o: ialglib.cpp ialglib.h
	g++  -c -O3 -g  -Wall -m64 ialglib.cpp
EntryCNV.o: EntryCNV.cpp EntryCNV.h
	g++  -c -O3 -g  -Wall -m64 EntryCNV.cpp
Reader.o: Reader.cpp Reader.h
	g++  -c  -O3 -g  -Wall -m64 Reader.cpp

Profiler.o: Profiler.cpp Profiler.h
	g++  -c  -O3 -g   -Wall -m64 Profiler.cpp
	
utils.o: utils.cpp utils.h
	g++  -c  -O3 -g   -Wall -m64 utils.cpp
	
segmentation.o: segmentation.cpp segmentation.h
	g++  -c  -O3 -g -Wall   -m64 segmentation.cpp
freec.o: freec.cpp freec.h
	g++  -c -O3 -g -Wall  -m64 freec.cpp
poissondistr.o: poissondistr.cpp poissondistr.h
	g++  -c -O3 -g -Wall  -m64 poissondistr.cpp
	
HMCan.o: HMCan.cpp HMCan.h
	g++  -c -O3 -g -Wall  -m64 HMCan.cpp

Parser.o: Parser.cpp Parser.h
	g++  -c -O3 -g -Wall  -m64 Parser.cpp

clean:
	rm -f HMCan *.o *~ *#












          
