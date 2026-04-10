# ChiCMaxima

![alt text](https://github.com/yousra291987/ChiCMaxima/blob/master/logo/logo-page-001.jpg)





ChiCMaxima: a pipeline for analyzing  and identificantion of chromation loops in CHi-C promoters data.

# Current release: ChiCMaxima 1.0

ChiCMaxima consists of 5 Rscripts: ChiCMaxima_Caller.r, ChiCMaxima_Collate.r, ChiCMaxima_MergeRep2.r,ChiCMaxima_MergeRepMany.r and ChiCMaxima_RepAnalysis.r, in addition to an R Browser (ChiCBrowser.r) to visualise CHi-C promoters data with called peaks. 



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

3/ ChiCMaxima_MergeRep2

Brief: Filters the called interactions from two biological replicates to only include those within a threshold genomic distance.
```
Usage (from folder containing scripts): Rscript ChiCMaxima_MergeRep2.r -a/--onepeak [INTERACTIONS FILE 1] -b/--twopeak [INTERACTIONS FILE 2] -d/--repdist [THRESHOLD DISTANCE] -o/--output [OUTPUT IBED FILE]

Arguments:

--onepeak, and
--twopeak. The output interaction files from ChiCMaxima_Caller, in the headed 12-column ibed format.
--repdist. The maximum genomic distance between the closest interactions within the two replicates for the interaction to be
maintained. Default: 0
--output. The file name for the output ibed file of the merged, filtered interactions list.
```


Outputs:

The merged ibed file (see above for column headings).

Example:
```R
Rscript ChiCMaxima_MergeRep2.r -a testdata/mESrep1_interactions.ibed -b testdata/mESrep2_interactions.ibed -d 20000 -o testdata/mEScombined_interactions.ibed
```
Returns 803 merged interactions.


4/ ChiCMaxima_MergeRepMany

Brief: Filters the called interactions from three or more biological replicates to only include those within a threshold distance.
```
Usage (from folder containing scripts): Rscript ChiCMaxima_MergeRepMany.r [THRESHOLD DISTANCE] [OUTPUT IBED FILE] [INTERACTIONS FILE 1] [INTERACTIONS FILE 2] [INTERACTIONS FILE 3] … [INTERACTIONS FILE n]

Arguments (note that these do not come with a “-“ or “--“ prefix):

THRESHOLD DISTANCE. As --repdist in ChiCMaxima_MergeRep2. There is no default setting, since all arguments must be entered.

OUTPUT IBED FILE. As --output in ChiCMaxima_MergeRep2.

INTERACTIONS FILES. A list of three or more interaction ibed files for merging.
```

Outputs:

The merged ibed file (see above for column headings).

Example:
```R
Rscript ChiCMaxima_MergeRepMany.r 20000 testdata/3combinedinteractions.ibed testdata/mESrep1_interactions.ibed testdata/mESrep2_interactions.ibed testdata/otherinteractions_chr15.ibed
```
Returns 643 merged interactions.

5/ ChiCMaxima_Collate

Brief: Reads in the input ibed files from different CHi-C experiments and collates them into one large ibed file with separate N columns for each dataset, suitable as input for ChiCBrowser.

```
Usage (from folder containing scripts): Rscript ChiCMaxima_Collate.r -k/--key [COLLATE KEY] -o/--output [OUTPUT IBED FILE]

Arguments:

--key. File name of a user-provided table which gives the names of all the input ibed files to be collated, along with their unique
identifiers for the experiment-specific N column. As shown in testdata/collate_key.txt, this is a non-headed 2-column table: IBED FILE,
IDENTIFIER.

--output. File name for the output collated ibed file.
```
Outputs:

Headed columns: ID_Bait, chr_Bait, start_Bait, end_Bait, Bait_name, ID_OE, chr_OE, start_OE, end_OE, OE_name, N.[IDENTIFIER 1], N.[IDENTIFIER 2], etc.

Example

```R
Rscript ChiCMaxima_Collate.r -k testdata/collate_key.txt -o testdata/testcollate.ibed
```
Outputs an ibed file, with N columns N.1 and N.2, identical to that already provided: testdata/mEScollated.ibed


# ChiCMaxima Browser

To run the browser, go to the folder containing the scripts and run the R environment. Then use the command: source(“ChiCBrowser.r”).
See AdditionalFile 2 for more explanation

![alt text](https://github.com/yousra291987/ChiCMaxima/blob/master/logo/screenshot_browser.png)


# ChiCMaxima Interactive Browser (R Shiny)

`app.R` provides an interactive web-based browser for visualising Promoter Capture Hi-C data, built with R Shiny. It reproduces all functionality of the original Tcl/Tk desktop GUI (`ChiCBrowser.r`) in a browser-accessible interface with no desktop windowing system required.

---

## Requirements

```r
install.packages(c("shiny", "data.table", "caTools"))
# Bioconductor packages
BiocManager::install(c("limma", "rtracklayer"))
```

---

## Usage

```r
shiny::runApp("app.R")
```

The app opens in your default browser. Follow the steps below.

---

## Workflow

### 1. Load interaction data
Upload an `.ibed` file via **"ibed interaction file"**. The file must contain the 10 standard ibed columns followed by one column per sample (replicate).

### 2. Assign replicates to conditions
After loading, a row appears for each sample column. Enter a positive integer to assign it to a condition (`1`, `2`, `3`, …) or `0` to exclude it. Click **Apply Conditions**.

### 3. Name and colour conditions
Each condition gets a default name (taken from the column header) and colour. Edit these under **Condition names & colours** and click **Apply Names & Colours**.

### 4. Navigate to a bait
Enter a bait gene name (must match `Bait_name` in the ibed file) and set either:
- a **window size** (bp centred on the bait fragment), or
- explicit **start / end coordinates**.

Click **Plot**.

### 5. Optional tracks

| Track type | Format | Notes |
|---|---|---|
| Gene annotation | `.txt` (tab-delimited, columns: `Chr`, `Start`, `End`, `Strand`, `Name`) | Drawn with directional arrows |
| Epigenomic signal | `.bw` (BigWig) or `.bedGraph` | Assign a **display level** (integer ≥ 1); tracks at the same level share a panel. Set to `0` to hide. |
| Called interactions | `.ibed` | Overlaid as points on the CHiC profile. Can be toggled on/off individually. |

### 6. Export
Click **Save Figure (PDF)** to download a publication-ready PDF of the current view.

---

## Plot description

The main panel shows a running-mean quantile-normalised CHiC interaction profile. Each condition is plotted as a separate coloured line. A dashed vertical line marks the bait fragment centre. Optional sub-panels below show (in order): called interaction overlays, gene annotation, and epigenomic signal tracks.

---

## Input file formats

**ibed file** – tab-delimited with a header row. The first 10 columns are standard ibed fields; subsequent columns are per-sample read counts:

```
chr_Bait  start_Bait  end_Bait  Bait_ID  Bait_name  chr_OE  start_OE  end_OE  OE_ID  OE_name  sample1  sample2  ...
```

**Gene annotation file** – tab-delimited with a header row:

```
Chr  Start  End  Strand  Name
```


# Dependencies:

    R version >= 3.2
    Bioconductor packages: GenomicRanges, limma, MASS, caTools, data.table, zoo, rtracklayer, psych, tcltk2, tkrplot 
 
 
 # Questions and contacts
 
 For FAQs, or for asking new questions, please contact: 
 
 yousra.ben-zouari@fmi.ch
 sexton@igbmc.fr

