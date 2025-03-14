# switch_error_screen

Haplotype-phased genomes can occasionally contain switch errors: an artificial switch in sequence of one haplotype to another within a contig that should only represent one haplotype. Switch errors can be detected by realigning the reads used for assembly to a phased contig or phased reference genome. Realignment of reads to a switch error junction should result in soft clipping of aligned reads on one side of the switch error where the change in contig sequence lacks read support.  

*switch_error_screen* is a bash script that uses soft-clipping information to flag potential switch errors in phased genomes and pangenomes. It can be applied to .BAM realignments of reads to a phased genome.

## Authors 
Samuel N. Bogan | University of California, Santa Cruz

Owen W. Moosman | University of California, Santa Cruz

Joanna L. Kelley | University of California, Santa Cruz

## Method

*switch_error_screen* reads in a BAM alignment of all hapmers (reads phased to a given haplotype) and unphased reads against a phased FASTA. Over a 10 kb sliding window, the script calculates an index of switch error by measuring three variables: the proportion of soft clipped reads within the window ($C$), the skewness of soft clipping toward 5’ and 3’ ends among these reads ($S$), and polarization of soft-clipping skewness ($P$). Each of these metrics are extracted from soft clipping data retained within the BAM CIGAR string. 

C is calculated as the number of soft clipped bases divided by the sum of all matched, inserted, and deleted bases. 

S is calculated as the variance in the per-read difference between the number of left and right soft-clipped bases. 

P is determined by the absolute difference between the number of reads clipped on the left versus the right, normalized by the total number of clipped reads.
 
Under switch error caused by incorrect phasing or misassembly, we expect to observe junctions at which raw reads exhibit (i) a high proportion of soft clipping ($C$), (ii) soft clipping on one side of each affected read ($S$), and (iii) soft clipping on the same side of a switch error junction among affected reads ($P$). The script then calculates the index of switch error risk according to the equation below: 
 
$$
Index = C(1+P)^2S
$$

Here, $(1+P)^2$ ensures that the index increases more rapidly as polarization deviates from zero. When polarization is zero ($P$=0), the term evaluates to 1 rather than 0, and does not eliminate risk associated with skewed soft-clipping being evident within a region.


## Usage

The script switch_error_screen.sh reads in two arguments: 

1. A. BAM file containing reads used for assembly, realigned to a phased assembly 

2. A .fai index of the phased assembly FASTA file

It can be run as follows:

    ./switch_error_screen.sh [bam] [fai]
    
The script will output the file 'soft_clip_metrics', which contain the soft clipping metrics described above over 10 kb sliding windows. 
    
## Citation

A manuscript reporting *switch_error_screen* is in review. For now, please cite this Github page and the listed authors if you use *switch_error_screen*.
    

