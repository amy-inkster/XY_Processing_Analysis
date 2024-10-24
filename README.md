# Who’s Afraid of the X? 
### Incorporating the X and Y chromosomes into the analysis of DNA methylation microarray data
This repository contains code from [Inkster et al. 2023](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-022-00477-0), with exemplars for processing and analysis (differential DNAme testing) available in the file ```aCA_XY_Processing_Analysis.md```.

Additional annotation files for the Illumina DNAme arrays indexing which CpG probes overlap repetitive elements, the X-transposed region, and cancer testis gene family members are also provided for both the 450K and EPIC arrays:
  1. Annotation_XY_450.csv
  2. Annotation_XY_850.csv
  
  To download these annotations, click on the file and  choose "Download" in the top right-hand corner, OR in R use:
  
  ```
  git450loc <- https://github.com/amy-inkster/XY_Processing_Analysis/blob/main/Annotation_XY_450.csv
  anno450 <- read.csv(git450loc)
  ```
  
  Please direct any feedback or questions to Amy Inkster at amymichelleinkster (at) gmail.com
