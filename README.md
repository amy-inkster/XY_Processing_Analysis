# XY_Processing_Analysis
This repository contains code from Inkster et al. 2022, with examplars for processing and analysis (differential DNAme testing) available in the file ```aCA_XY_Processing_Analysis.md```.

Additional annotation files for the Illumina DNAme arrays indexing which CpG probes overlap repetitive elements, the X-transposed region, and cancer testis gene family members are also provided for both the 450K and EPIC arrays:
  1. Annotation_XY_450.csv
  2. Annotation_XY_850.csv
  
  To download these annotations, click on the file and  choose "Download" in the top right-hand corner, OR in R use:
  
  ```
  git450loc <- https://github.com/amy-inkster/XY_Processing_Analysis/blob/main/Annotation_XY_450.csv
  anno450 <- read.csv(git450loc)
  ```
