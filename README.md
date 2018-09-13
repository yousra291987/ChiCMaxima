# ChiCMaxima

ChiCMaxima pipeline for analyzing  and identificantion of chromation loops in CHi-C promoters data. ChiCMaxima consists of 5 Rscripts: ChiCMaxima_Caller.r, ChiCMaxima_Collate.r, ChiCMaxima_MergeRep2.r,ChiCMaxima_MergeRepMany.r and ChiCMaxima_RepAnalysis.r, in addition to an R Browser (ChiCBrowser.r) to visualise CHi-C promoters data with called peaks. 

Current release: ChiCMaxima 1.0

# Input format

The current version uses ibed matrices with 11 columns in which the coordinates of both interacting regions plus the number of reads are present.

 For example:
 
|ID_Bait|chr_Bait|start_Bait|end_Bait|Bait_name|ID_OE|chr_OE|start_OE|end_OE|OE_name|N|
|-------|--------|----------|--------|---------|-----|------|--------|------|-------|-|
|30|chr1|3090913|3092556|U6.149|590|chr1|4592259|4592779|.|0|
|30|chr1|3090913|3092556|U6.149|593|chr1|4595997|4596467|.|1|
|30|chr1|3090913|3092556|U6.149|596|chr1|4605050|4610398|.|2|


Scripts (some from CHiCAGO) for generating and reconverting these files…


# Usage

1/ ChiCMaxima_Caller

Brief: Calls CHi-C interactions from single datasets, based on local maximum computation and approximate prediction of cis-decay.

```
Usage (from folder containing scripts): Rscript ChiCMaxima_Caller.r -i/--input [INPUT IBED FILE] -o/--output [OUTPUT PREFIX]
-w/--window_size [LOCAL MAXIMUM CALLING WINDOW] -s/--loess_span [LOESS SPAN] -c/--cis_window [GENOMIC SEPARATION THRESHOLD]

Arguments:
 --input. 11-column ibed input file (headed or not; see section 2 for format).
 --output. Folder and prefix for output files (called interactions, and list of non-assessed baits).
 --window_size. Number of covered restriction fragments within which the local maximum is called in sliding windows. 
                 Default:50
 --loess_span. The loess span parameter for smoothing the virtual 4C profile. Default: 0.05
 --cis_window. Interactions are assessed within this many bp of the bait start coordinate. Default: 1500000
```



Outputs:

  . output_interactions.ibed. The called interactions, in the same ibed format as the input, with a twelfth column, Enrichment - log2(local maximum score/fitted cis-decay score).

  . output_poorlycovered.txt.A list of the baits with sufficient coverage for local maximum assessment. This and their numbers can be used to inform the user of which -w parameter to use.

Example:
```R
Rscript ChiCMaxima_Caller.r -i testdata/mESrep1_chr15.ibed -o testdata/mESrep1 -w 20
```
Returns 2561 interactions, with 33 poorly covered baits.
```R
Rscript ChiCMaxima_Caller.r -i testdata/mESrep2_chr15.ibed -o testdata/mESrep2 -w 20
```
Returns 1645 interactions, with 65 poorly covered baits.

2/ ChiCMaxima_RepAnalysis

Brief: Computes the distributions of closest distances between called interactions within two outputs of CHiCMaxima_Caller. The goal is to help the user choose an appropriate -d setting for ChiCMaxima_MergeRep2/ChiCMaxima_MergeRepMany (see section 5).
```
Usage (from folder containing scripts): Rscript ChiCMaxima_RepAnalysis.r -a/--file1 [INTERACTIONS FILE 1] -b/--file2 [INTERACTIONS FILE 2] -o/--output [OUTPUT PREFIX]

Arguments:

--file1, and

--file2. The output interaction files from ChiCMaxima_Caller, in the headed 12-column ibed format.

--output. Folder and prefix for output files (histogram, cumulative frequency plot, quantile table).
```

Outputs:

  . output_hist.png. The histogram for the closest distance distribution.

  . output_ecdf.png. The cumulative frequency plot for closest distances.

  . output_percentiles.txt. The table giving the genomic distance for every fifth percentile.

Example:
```R
Rscript ChiCMaxima_RepAnalysis.r -a testdata/mESrep1_interactions.ibed -b testdata/mESrep2_interactions.ibed -o testdata/mEScombined.
```
Returns the histogram, cumulative frequency plot and percentiles table (25th percentile is 0; median is 18687 bp; 75th percentile is 79633.5 bp).








































# ChiCMaxima Browser

Example of CHi-C promoters data visualized by the ChiCMaxima Browser

# Dependencies:

    R version >= 3.2
    Bioconductor packages: GenomicRanges, limma, MASS, caTools, data.table, zoo, rtracklayer, psych, tcltk2, tkrplot 
 
 
 # Questions and contacts
 
 For FAQs, or for asking new questions, please see our forum: 
 
 yousra.ben-zouari@fmi.ch
 sexton@igbmc.fr

