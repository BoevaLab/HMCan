gccount: main.o ConfigFile.o Chameleon.o myFunc.o ChrCopyNumber.o GenomeCopyNumber.o
	g++ -O3 -Wall -m64 -o gccount main.o ConfigFile.o Chameleon.o myFunc.o ChrCopyNumber.o GenomeCopyNumber.o
	rm main.o
	rm ConfigFile.o
	rm Chameleon.o	
	rm myFunc.o	
	rm ChrCopyNumber.o 
	rm GenomeCopyNumber.o		
main.o: main.cpp SVfinder.h
	g++ -g -O3 -c -Wall -m64 main.cpp
ConfigFile.o: ConfigFile.cpp ConfigFile.h
	g++ -g -O3 -c -Wall -m64  ConfigFile.cpp
Chameleon.o: Chameleon.cpp Chameleon.h
	g++ -g -O3 -c -Wall -m64  Chameleon.cpp
myFunc.o: myFunc.cpp myFunc.h
	g++ -g -c -O3 -Wall -m64  myFunc.cpp	
GenomeCopyNumber.o: GenomeCopyNumber.cpp GenomeCopyNumber.h
	g++ -g -c -O3 -Wall -m64  GenomeCopyNumber.cpp	
ChrCopyNumber.o: ChrCopyNumber.cpp ChrCopyNumber.h
	g++ -g -c -O3 -Wall -m64  ChrCopyNumber.cpp	
clean:
	rm -f gccount *.o *~ *#          
