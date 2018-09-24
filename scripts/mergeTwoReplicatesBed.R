library(ucscGenomes())

args <- commandArgs()

if (length(args)<6) {
  print ("cat mergeTwoReplicatesBed.R | R --slave --args < bed1 > < bed2 > < outfile.bed > [<minimal score for short peaks>]")
  quit()
}

bedfile1=args[4]
bedfile2=args[5]
outfile=args[6]
minscore=-100;

if (length(args)==7) {
  minscore= type.convert(args[8])
}

final.peak1=import.bed(bedfile1)
final.peak2=import.bed(bedfile2)
common.peaks <- intersect(final.peak1,final.peak2) 
all.peaks <- append(final.peak1,final.peak2) 
all.peaks <- reduce(all.peaks)

#get close overlapping peaks:
mergedSet <- common.peaks[1]
for (i in c(1:(length(common.peaks)-1))) {
  a <- common.peaks[i]
  b <- common.peaks[i+1]
  overlap <- subsetByOverlaps(all.peaks,c(a,b)) 
  if (length(overlap)>1) { #save two peaks separately
    mergedSet <- append(mergedSet,append(a,b))   
  } else {
    mergedSet <- append(mergedSet,GRanges(chrom(a),IRanges(start(a),end(b))))
  }
  if (i%%1000==0) {print(paste("processed" ,i, "peaks", sep=" "))}
}
mergedSet <- reduce(mergedSet)

highPeaks.1 <- subsetByOverlaps(final.peak1,mergedSet )[which(score(subsetByOverlaps(final.peak1,mergedSet )) > minscore)]
highPeaks.2 <- subsetByOverlaps(final.peak2,mergedSet )[which(score(subsetByOverlaps(final.peak2,mergedSet )) > minscore)]
highPeaks.1.2 <- append(highPeaks.1,highPeaks.2)
high.mergedSet <- subsetByOverlaps(mergedSet,highPeaks.1.2)

export(high.mergedSet,BEDFile(outfile), format = "bed")