# ChiCMaxima

ChiCMaxima pipeline for analyzing  and identificantion of chromation loops in CHi-C promoters data.
Current release: ChiCMaxima 0.9


# Installation


# Input format

The current version uses ibed matrices with 11 columns in which the coordinates of both interacting regions plus the number of reads are present.
 For example:
 
|ID_Bait|chr_Bait|start_Bait|end_Bait|Bait_name|ID_OE|chr_OE|start_OE|end_OE|OE_name|N|
|-------|--------|----------|--------|---------|-----|------|--------|------|-------|-|
|30|chr1|3090913|3092556|U6.149|590|chr1|4592259|4592779|.|0|
|30|chr1|3090913|3092556|U6.149|593|chr1|4595997|4596467|.|1|
|30|chr1|3090913|3092556|U6.149|596|chr1|4605050|4610398|.|2|

ChiCMaxima will produce files with 11 columns that look like the Input which corresponds to the significant interactions.

# Usage

You can get a list of ChiCMaxima' command line arguments by passing the --help parameter. The current arguments are:

 Options:
    
    -i/--input			     <string>	   [default:input  ]
    -o/--output       <string>    [default:./ouput]
    -w/--window_size  <string>    [default: 50    ]
    -s/--loess_span   <string>    [default: 0.05  ]
    -h/--help                       
    

An example run on chromosome 1 of a mouse embryonic stem cells (Schoenfelder et al, 2015):

ChiCMaxima_caller.r -i Chr1_ES.ibed -o Chr1_ES_significant_interactions.ibed


We have also added a simple example in the "examples/" directory for your convenience.

# Output format

ChiCMaxima will produce files containing the significant interactions with the same columns as the input table.

# Dependencies:

    R version >= 3.2
    Bioconductor packages: Rsamtools,GenomicRanges, limma, MASS, caTools, data.table, base, zoo, RcppRoll, psych, plyr, ICSNP 

