##HMCan: detection of histone modifications in cancer samples

**Authors:** Haitham Ashoor, Valentina Boeva
**Contact:** haitham.ashoor @ jax . org, Valentina.Boeva @ inserm . fr

###Compilation:
In order to compile HMCan on your system use the provided make file with HMCan
In order to complie the HMCan, write the following command:

`make`

###Running HMCan:
You can run HMCan using the following command:

`./HMCan <TargetFile> <ControlFile> <Configuration file> <Name>`

If you do not have Input data, or your input file is a .WIG, use can run HMCan without a control file:


`./HMCan <TargetFile> - <Configuration file> <Name>`
*Do not forget the dash!*

**TargetFile**: Alignment file for target ChIP-seq data.
**ControlFile**: Alignment file for control ChIP-seq data.
**Configuration file**: Configuration file. See mandatory and optional parameters below.
**Name**: prefix name for the output files. All output files will have this name as prefix. Add here the output directory path as well!


####Configuration file:
In configuration file you can provide parameters for HMCan algorithm,
the description of the parameters is described below.

***Default parameters are provided in /configurations/HMCan.config.*.txt**

**Parameters description:**

Parameter | Mandatory/Optional | Description
--------------|---------------------------|----------------
format | Optional. Default: BAM | HMCan accepts BAM, SAM and BED alignment formats, or WIG files with pre-calculated and normalized densities in the "variableStep" or "fixedStep" format. **In order to process BAM file, `samtools` should be installed on the system.**
GCIndex | Mandatory | This file contains GC content and mapability scores for large regions of DNA. The file should be formated in .cnp format. We provide such files in **/data/** for mm9, mm10, hg19 and hg38 assemblies. Please set **largeBinLength** to the value indicated in the corresponding GC_profile file (25000 or 100000).
smallBinLength | Mandatory | Bin size to be used to construct signal profiles. Use 50 for ChIP-seq of Histone marks, and 10 for TFs or ATAC-seq. *When using .wig as input, provide the correct step length here*
largeBinLength | Mandatory | Bin size to be used to calculate copy number. Please note that should be the same length for GCIndex
genomePath | Mandatory | A path to chromosome sequences in .fa format. Each chromosome should be located in a separate file. One can download hg19 chromosomes from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/
minLength | Optional. Default: 145. **Use only for ChIP-seq Single End data** | Minimum fragment length used in the ChIP-seq experiment
medLength | Optional. Default: 150. **Use only for ChIP-seq Single End data** | Median fragment length used in the ChIP-seq experiment
maxLength | Optional. Default: 155. **Use only for ChIP-seq Single End data** | Maximum fragment length used in the ChIP-seq experiment
pvalueThreshold | Optional. Default: 0.01 | p-value threshold of Poisson's single side exact test. Use 0.05 for noisy data to be more sensitive in peak detection.
mergeDistance | Optional. Default: 200 | max distance to merge single peaks into region
iterationThreshold | Optional. Default: 5 | Score threshold to remove peaks in the iteration stage
finalThreshold | Optional. Default: 0 | Score threshold to report peaks or regions
maxIter | Optional. Default: 10 | Maximum number of iterations for HMCan algorithm
printWig | Optional. Default: False | A boolean option enables the user to print density in WIG files or not. **Set to 'True' if you want to caclulate fragment density
printPosterior | Optional. Default: False | A boolean option enables the user to print the bins posterior probabilities in a WIG format
PorsteriorProb | Optional. Default: 0.7 | Threshold for posterior probability to consider bin to be in the "peak" state
blackListFile | Optional | A file that contains regions to be ignored by HMCan. It should be a three column file chr,start,end. Example for Human blacklist file is in data folder
printBEDGraph | Optional. Default: False | A boolean option enables the user to print bedgraph files for normalized profiles of ChIP and Input libraries separately
CallPeaks | Optional. Default: True | A boolean option enables the user either to perform peak calling task of HMCan or just print normalized profiles only
pairedEnds | Optional. Default: False | Boolean. Set to 'True' when you have paired-end data. In this case, you can ignore parameters 4-6: minLength, medLength & maxLength
Type | Optional. Default: "Type ChIP-seq". | Set this parameter to 'ATAC-seq' if you are analyzing ATAC-seq data. ATAC-seq data will not be normalized for the GC-content.
GCmergeDistance | Optional. Default: 1000 | Distance (in bp) to remove closely located peaks during the calculation of GC-content bias in ChIP data
RemoveDuplicates | Optional. Default: True. | Set to 'False' if you do not want to remove duplicate reads. Can be important to detect signal in the amplification regions (copy number >50).
**calculateEmpiricalPvalue** | Optional. Default: False | Set to 'True' to caclulate empirical p-value using peaks in the Input. Providing the input dataset is mandatory to use this option.
ifKeepChIPdensityIntact | Optional. Default: True | Set to 'False' if you want to get HMCan-v.1.39-like behaviour where ChIP-seq density can be artificially increased when there are more reads in the Input .BAM file compared to ChIP.
CNAnormalization | Optional. Default: True | Set to 'False' if you do not want to normalize the ChIP-seq/ATAC-seq signal for copy number alterations, e.g., when your sample is diploid and not coming from a tumor or cancer cell line.
genomeLengthFile | Mandatory when reading .WIG | Provide a .fai file (e.g. [link](https://github.com/reineckef/quandico/blob/master/RPackage/hg19.fa.fai))

####Output files:
HMCan produces three main output files for narrow **peaks**, larger **regions** with the signal, and whole genome density.

Peaks file is in .narrowPeak format (https://genome.ucsc.edu/FAQ/FAQformat.html#format12) and regions file is in the .BED format.
Density files are in wig format for final normalized density and BedGraph for ChIP and Input normalized densities.
You can use wigToBigWig Utility to transform .wig to .bw: [link](http://genome.ucsc.edu/goldenPath/help/bigWig.html)

Output name | Optional or default | Format | Comment
------------------|---------------|-------------|------------
**_regions.bed**| Default | Chromosome, Start, End, PeakOrder,  Score, strand (.), empirical p-value | Large regions with signal in .BED format. p-value is provided if **'calculateEmpiricalPvalue True'**
**_peaks.narrowPeak**| Default | Chromosome, Start, End, PeakOrder,  Score, strand (.), peak hight, -log10(empirical p-value), -1, relative position of peak maximum | Narrow regions with signal in .narrowPeak format. p-value is provided if **'calculateEmpiricalPvalue True'**
**.wig** | Optional. Set 'printWig True' | wig | Normalized density profile
**_CNV_profile.txt** | Default | ControlFREEC _ratio.txt format | Copy number profile from Input or ChIP (if Input is not profided). Can be visualized with http://boevalab.com/FREEC/tutorial.html#VISUAL


####Citation:
If you use HMCan please cite "HMCan: a method for detecting chromatin modifications in cancer samples using ChIP-seq data"
Haitham Ashoor; Aurelie Herault; Aurelie Kamoun; Francois Radvanyi; Vladimir B. Bajic; Emmanuel Barillot; Valentina Boeva.
Bioinformatics, 2013, 29(23):2979-86. doi: 10.1093/bioinformatics/btt524

**NOTE:** starting from V1.16 chromosome names in the GC index file should match exactly names in alignment and provided genome file.
