#!/bin/bash

# Check for input arguments
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <input.bam> <reference.fai>"
  exit 1
fi

# Input BAM file and reference FAI file
BAM=$1
FAI=$2

# Output files
OUTPUT="soft_clip_metrics.tsv"
ERROR_FILE="errors.log"

# Define sliding window size (1 kb)
WINDOW_SIZE=10000

# Initialize error file and output file
echo "Starting script at $(date)" > "$ERROR_FILE"
echo -e "chromosome\tstart\tend\tP_s\tS\tR\tW" > "$OUTPUT"

# Get chromosomes present in the BAM file
CHROMOSOMES=$(samtools idxstats "$BAM" 2>>"$ERROR_FILE" | awk '$3 > 0 {print $1}')

# Get the total number of windows for the loading bar
TOTAL_WINDOWS=0
for chrom in $CHROMOSOMES; do
  length=$(grep -w "$chrom" "$FAI" 2>>"$ERROR_FILE" | cut -f2)
  if [ -n "$length" ]; then
    TOTAL_WINDOWS=$((TOTAL_WINDOWS + (length + WINDOW_SIZE - 1) / WINDOW_SIZE))
  fi
done

# Progress bar variables
PROGRESS=0
update_progress_bar() {
  local percent=$((100 * PROGRESS / TOTAL_WINDOWS))
  echo -ne "\rProgress: ["
  for ((i = 0; i < percent / 2; i++)); do echo -n "#"; done
  for ((i = percent / 2; i < 50; i++)); do echo -n " "; done
  echo -n "] $percent%"
}

# Loop through each chromosome present in the BAM file
for chrom in $CHROMOSOMES; do
  length=$(grep -w "$chrom" "$FAI" 2>>"$ERROR_FILE" | cut -f2)
  if [ -z "$length" ]; then
    echo "Warning: Chromosome $chrom not found in reference FAI file. Skipping." >>"$ERROR_FILE"
    continue
  fi

  # Loop through windows in the chromosome
  for ((start = 0; start < length; start += WINDOW_SIZE)); do
    end=$((start + WINDOW_SIZE - 1))
    if [ "$end" -ge "$length" ]; then
      end=$((length - 1))
    fi

    # Extract alignments within the window
    samtools view -F 4 "$BAM" "$chrom:$start-$end" 2>>"$ERROR_FILE" >alignments.tmp

    # Calculate soft-clipped bases (P_s)
    total_bases=0
    soft_clipped_bases=0
    awk '{
      cigar=$6;
      while (cigar ~ /[0-9]+[SMID]/) {
        match(cigar, /[0-9]+[SMID]/);
        len = substr(cigar, RSTART, RLENGTH - 1);
        op = substr(cigar, RSTART + RLENGTH - 1, 1);
        cigar = substr(cigar, RSTART + RLENGTH);
        if (op == "S") {
          soft_clipped_bases += len;
        }
        if (op == "M" || op == "I" || op == "D") {
          total_bases += len;
        }
      }
    } END {
      if (total_bases > 0) {
        print soft_clipped_bases / total_bases;
      } else {
        print 0;
      }
    }' alignments.tmp >Ps.tmp 2>>"$ERROR_FILE"
    P_s=$(cat Ps.tmp)

    # Calculate skewness (S) of soft clipping
    awk '{
      cigar=$6;
      left=0;
      right=0;
      if (cigar ~ /^[0-9]+S/) {
        match(cigar, /^[0-9]+S/);
        left = substr(cigar, RSTART, RLENGTH - 1);
      }
      if (cigar ~ /[0-9]+S$/) {
        match(cigar, /[0-9]+S$/);
        right = substr(cigar, RSTART, RLENGTH - 1);
      }
      print left - right;
    }' alignments.tmp >skewness.tmp 2>>"$ERROR_FILE"
    S=$(awk '{sum+=$1; sumsq+=$1*$1} END {print (NR>0 ? (sumsq/NR - (sum/NR)^2) : 0)}' skewness.tmp 2>>"$ERROR_FILE")

    # Calculate polarization (R)
    R=$(awk '{
      if ($1 > 0) left++;
      else if ($1 < 0) right++;
    } END {
      if (left + right > 0) {
        diff = left - right;
        print (diff < 0 ? -diff : diff) / (left + right); # Manual absolute value
      } else {
        print 0;
      }
    }' skewness.tmp 2>>"$ERROR_FILE")

    # Calculate width penalty (W)
    W=$((end - start + 1))

    # Append results to the output file
    echo -e "$chrom\t$start\t$end\t$P_s\t$S\t$R\t$W" >>"$OUTPUT"

    # Update progress
    PROGRESS=$((PROGRESS + 1))
    update_progress_bar
  done
done

# Cleanup temporary files
rm -f alignments.tmp Ps.tmp skewness.tmp

# Finalize progress bar
echo -e "\nDone! Results saved to $OUTPUT, errors logged in $ERROR_FILE"
