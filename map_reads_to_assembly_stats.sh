#!/bin/bash


# Usage: ./map_reads_to_assembly_stats.sh <input.bam> <output.tsv> [threads]

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <input.bam> <output.tsv> [threads]"
  exit 1
fi

BAM_FILE=$1
OUTPUT_TSV=$2
# Default to 4 threads if not specified
THREADS=${3:-4}

echo "Using $THREADS threads for processing (where supported)"

# Ensure the BAM file is sorted and indexed
echo "Ensuring BAM file is sorted and indexed..."
# sort does support threading
samtools sort -@ "$THREADS" -o "${BAM_FILE%.bam}.sorted.bam" "$BAM_FILE"
# index in v1.12 doesn't support threading
samtools index "${BAM_FILE%.bam}.sorted.bam"
SORTED_BAM="${BAM_FILE%.bam}.sorted.bam"

# Initialize TSV file with headers
echo -e "Metric\tValue" > "$OUTPUT_TSV"

# Overall mapping statistics
echo "Calculating overall mapping statistics..."
# view supports threading
TOTAL_READS=$(samtools view -@ "$THREADS" -c "$SORTED_BAM")
MAPPED_READS=$(samtools view -@ "$THREADS" -c -F 4 "$SORTED_BAM")

# Fix: Prevent division by zero
if [ "$TOTAL_READS" -gt 0 ]; then
  MAPPING_RATE=$(awk "BEGIN {printf \"%.2f\", ($MAPPED_READS/$TOTAL_READS)*100}")
else
  MAPPING_RATE="0.00"
fi

echo -e "Total_Reads\t$TOTAL_READS" >> "$OUTPUT_TSV"
echo -e "Mapped_Reads\t$MAPPED_READS" >> "$OUTPUT_TSV"
echo -e "Mapping_Rate_Percent\t$MAPPING_RATE" >> "$OUTPUT_TSV"

# Per-reference statistics
echo "Calculating per-reference statistics..."
echo -e "\n# Per-Reference Statistics" >> "$OUTPUT_TSV"
echo -e "Reference\tLength\tMapped_Reads\tCoverage\tZero_Cov_Bases\tZero_Cov_Percent" >> "$OUTPUT_TSV"

# Create temp file for coverage stats
COVERAGE_TEMP=$(mktemp)
# coverage in v1.12 might not support threading
samtools coverage "$SORTED_BAM" > "$COVERAGE_TEMP"

# Process each reference
tail -n +2 "$COVERAGE_TEMP" | while read -r line; do
  REF=$(echo "$line" | awk '{print $1}')
  START=$(echo "$line" | awk '{print $2}')
  END=$(echo "$line" | awk '{print $3}')
  LENGTH=$((END - START + 1))
  MAPPED=$(echo "$line" | awk '{print $4}')
  COV=$(echo "$line" | awk '{print $7}')
  
  # Calculate zero coverage positions
  DEPTH_TEMP=$(mktemp)
  # depth in v1.12 might not support threading
  samtools depth -a -r "$REF" "$SORTED_BAM" > "$DEPTH_TEMP"
  
  # Fix: Properly count lines with zero coverage
  # The format of samtools depth output is: <ref> <pos> <depth>
  ZERO_COV=$(awk '$3 == 0 {count++} END {print count}' "$DEPTH_TEMP")
  
  # Fix: Handle case when ZERO_COV is empty
  if [ -z "$ZERO_COV" ]; then
    ZERO_COV=0
  fi
  
  # Fix: Prevent division by zero
  if [ "$LENGTH" -gt 0 ]; then
    ZERO_PERCENT=$(awk "BEGIN {printf \"%.2f\", ($ZERO_COV/$LENGTH)*100}")
  else
    ZERO_PERCENT="0.00"
  fi
  
  echo -e "$REF\t$LENGTH\t$MAPPED\t$COV\t$ZERO_COV\t$ZERO_PERCENT" >> "$OUTPUT_TSV"
  rm "$DEPTH_TEMP"
done

# Identify zero coverage regions (gaps) using samtools only
echo "Finding zero coverage regions..."
echo -e "\n# Zero Coverage Regions (>20bp)" >> "$OUTPUT_TSV"
echo -e "Reference\tStart\tEnd\tLength" >> "$OUTPUT_TSV"

# Get sequence lengths
SEQLEN_TEMP=$(mktemp)
samtools view -H "$SORTED_BAM" | grep "@SQ" | sed 's/.*SN:\([^\t]*\)\tLN:\([0-9]*\).*/\1\t\2/' > "$SEQLEN_TEMP"

# Process each reference sequence for gaps
while read -r ref len; do
  # Use samtools depth with -a to report all positions, including zero coverage
  samtools depth -a -r "$ref" "$SORTED_BAM" | 
    awk -v ref="$ref" -v output="$OUTPUT_TSV" '
    BEGIN {OFS="\t"; gap_start=0; in_gap=0}
    {
      if ($3 == 0) {
        # Zero coverage position
        if (!in_gap) {
          # Start of a new gap
          gap_start = $2
          in_gap = 1
        }
      } else {
        # Coverage > 0
        if (in_gap) {
          # End of a gap
          gap_end = $2 - 1
          gap_len = gap_end - gap_start + 1
          
          # Only report gaps > 20bp
          if (gap_len > 20) {
            print ref, gap_start, gap_end, gap_len >> output
          }
          in_gap = 0
        }
      }
      prev_pos = $2
    }
    END {
      if (in_gap) {
        gap_end = prev_pos
        gap_len = gap_end - gap_start + 1
        if (gap_len > 20) {
          print ref, gap_start, gap_end, gap_len >> output
        }
      }
    }'
done < "$SEQLEN_TEMP"

# Coverage distribution analysis
echo -e "\n# Coverage Distribution" >> "$OUTPUT_TSV"
echo -e "Coverage_Range\tBases\tPercent_of_Genome" >> "$OUTPUT_TSV"

# Get the total genome size
GENOME_SIZE=$(awk '{sum+=$2} END {print sum}' "$SEQLEN_TEMP")
DEPTH_ALL_TEMP=$(mktemp)

# Get depth for all positions
echo "Generating genome-wide coverage distribution..."
samtools depth -a "$SORTED_BAM" > "$DEPTH_ALL_TEMP"

# Define coverage ranges to analyze
COV_RANGES=("0" "1-10" "11-20" "21-50" "51-100" "101-200" "201+")

# Count bases in each coverage range
# Fix: Correctly count 0x coverage bases
COV_0=$(awk '$3 == 0 {count++} END {print count}' "$DEPTH_ALL_TEMP")

# Handle case where depth file might not include all zero coverage positions
if [ -z "$COV_0" ]; then
  COV_0=0
fi

# If COV_0 is still 0, calculate it as the difference between genome size and total reported positions
if [ "$COV_0" -eq 0 ]; then
  TOTAL_REPORTED=$(wc -l < "$DEPTH_ALL_TEMP")
  if [ "$TOTAL_REPORTED" -lt "$GENOME_SIZE" ]; then
    COV_0=$((GENOME_SIZE - TOTAL_REPORTED))
  fi
fi

COV_1_10=$(awk '$3>=1 && $3<=10' "$DEPTH_ALL_TEMP" | wc -l)
COV_11_20=$(awk '$3>=11 && $3<=20' "$DEPTH_ALL_TEMP" | wc -l)
COV_21_50=$(awk '$3>=21 && $3<=50' "$DEPTH_ALL_TEMP" | wc -l)
COV_51_100=$(awk '$3>=51 && $3<=100' "$DEPTH_ALL_TEMP" | wc -l)
COV_101_200=$(awk '$3>=101 && $3<=200' "$DEPTH_ALL_TEMP" | wc -l)
COV_201_PLUS=$(awk '$3>=201' "$DEPTH_ALL_TEMP" | wc -l)

# Calculate percentages - Fix: Prevent division by zero
if [ "$GENOME_SIZE" -gt 0 ]; then
  PCT_0=$(awk "BEGIN {printf \"%.2f\", ($COV_0/$GENOME_SIZE)*100}")
  PCT_1_10=$(awk "BEGIN {printf \"%.2f\", ($COV_1_10/$GENOME_SIZE)*100}")
  PCT_11_20=$(awk "BEGIN {printf \"%.2f\", ($COV_11_20/$GENOME_SIZE)*100}")
  PCT_21_50=$(awk "BEGIN {printf \"%.2f\", ($COV_21_50/$GENOME_SIZE)*100}")
  PCT_51_100=$(awk "BEGIN {printf \"%.2f\", ($COV_51_100/$GENOME_SIZE)*100}")
  PCT_101_200=$(awk "BEGIN {printf \"%.2f\", ($COV_101_200/$GENOME_SIZE)*100}")
  PCT_201_PLUS=$(awk "BEGIN {printf \"%.2f\", ($COV_201_PLUS/$GENOME_SIZE)*100}")
else
  PCT_0="0.00"
  PCT_1_10="0.00"
  PCT_11_20="0.00"
  PCT_21_50="0.00"
  PCT_51_100="0.00"
  PCT_101_200="0.00"
  PCT_201_PLUS="0.00"
fi

# Add to output
echo -e "0x\t$COV_0\t$PCT_0%" >> "$OUTPUT_TSV"
echo -e "1-10x\t$COV_1_10\t$PCT_1_10%" >> "$OUTPUT_TSV"
echo -e "11-20x\t$COV_11_20\t$PCT_11_20%" >> "$OUTPUT_TSV"
echo -e "21-50x\t$COV_21_50\t$PCT_21_50%" >> "$OUTPUT_TSV"
echo -e "51-100x\t$COV_51_100\t$PCT_51_100%" >> "$OUTPUT_TSV"
echo -e "101-200x\t$COV_101_200\t$PCT_101_200%" >> "$OUTPUT_TSV"
echo -e "201+x\t$COV_201_PLUS\t$PCT_201_PLUS%" >> "$OUTPUT_TSV"

# Calculate overall statistics
TOTAL_COVERED=$((GENOME_SIZE - COV_0))
if [ "$GENOME_SIZE" -gt 0 ]; then
  PCT_COVERED=$(awk "BEGIN {printf \"%.2f\", ($TOTAL_COVERED/$GENOME_SIZE)*100}")
else
  PCT_COVERED="0.00"
fi
echo -e "Total_Covered_Bases\t$TOTAL_COVERED\t$PCT_COVERED%" >> "$OUTPUT_TSV"

# Calculate mean coverage of covered bases (excluding zeros)
# Fix: Prevent division by zero in mean calculation
MEAN_COV_NONZERO=$(awk '$3>0 {sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print "0.00"}' "$DEPTH_ALL_TEMP")
echo -e "Mean_Coverage_NonZero\t$MEAN_COV_NONZERO\tN/A" >> "$OUTPUT_TSV"

# Add circular genome detection
echo -e "\n# Circularity Analysis" >> "$OUTPUT_TSV"
echo -e "Reference\tCircularity_Score\tStatus\tEvidence" >> "$OUTPUT_TSV"

echo "Checking for circular genome evidence..."

while read -r ref len; do
  echo "Analyzing circularity for $ref..."
  
  # Skip small contigs
  if [ "$len" -lt 20000 ]; then
    echo -e "$ref\t0\tToo_Small\tLength < 20kb" >> "$OUTPUT_TSV"
    continue
  fi
  
  # Define regions to check (first and last 1000bp)
  START_REGION="${ref}:1-1000"
  END_REGION="${ref}:$((len-1000+1))-$len"
  
  # Create temp files for analysis
  START_READS=$(mktemp)
  END_READS=$(mktemp)
  SPANNING_READS=$(mktemp)
  
  # Extract reads mapping to start and end regions
  # view does support threading in v1.12
  samtools view -@ "$THREADS" "$SORTED_BAM" "$START_REGION" | cut -f1 > "$START_READS"
  samtools view -@ "$THREADS" "$SORTED_BAM" "$END_REGION" | cut -f1 > "$END_READS"
  
  # Find reads that map to both start and end (potential spanning reads)
  sort "$START_READS" "$END_READS" | uniq -d > "$SPANNING_READS"
  SPAN_COUNT=$(wc -l < "$SPANNING_READS")
  
  # Get coverage at start and end compared to middle
  # depth doesn't support -@ in v1.12
  START_COV=$(samtools depth -r "$START_REGION" "$SORTED_BAM" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
  END_COV=$(samtools depth -r "$END_REGION" "$SORTED_BAM" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
  MIDDLE_REGION="${ref}:$((len/2-500)):$((len/2+500))"
  MIDDLE_COV=$(samtools depth -r "$MIDDLE_REGION" "$SORTED_BAM" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
  
  # Calculate circularity metrics (avoid division by zero)
  if [ "$MIDDLE_COV" != "0" ] && [ ! -z "$MIDDLE_COV" ]; then
    START_END_RATIO=$(awk "BEGIN {print (($START_COV+$END_COV)/2)/$MIDDLE_COV}")
  else
    START_END_RATIO=0
  fi
  
  # Look for supplementary alignments at start/end (split reads)
  # view supports threading
  SPLIT_COUNT=$(samtools view -@ "$THREADS" -f 2048 "$SORTED_BAM" "$START_REGION" "$END_REGION" | wc -l)
  
  # Calculate circularity score (0-100)
  # Factors: spanning reads, coverage consistency, split reads
  CIRC_SCORE=0
  EVIDENCE=""
  
  # Add score for spanning reads (up to 40 points)
  if [ "$SPAN_COUNT" -ge 10 ]; then
    CIRC_SCORE=$((CIRC_SCORE + 40))
    EVIDENCE="${EVIDENCE}Spanning_reads:${SPAN_COUNT};"
  elif [ "$SPAN_COUNT" -ge 5 ]; then
    CIRC_SCORE=$((CIRC_SCORE + 20))
    EVIDENCE="${EVIDENCE}Spanning_reads:${SPAN_COUNT};"
  elif [ "$SPAN_COUNT" -ge 1 ]; then
    CIRC_SCORE=$((CIRC_SCORE + 10))
    EVIDENCE="${EVIDENCE}Spanning_reads:${SPAN_COUNT};"
  fi
  
  # Add score for coverage consistency (up to 40 points)
  if (( $(awk "BEGIN {print ($START_END_RATIO > 0.9 && $START_END_RATIO < 1.1)}") )); then
    CIRC_SCORE=$((CIRC_SCORE + 40))
    EVIDENCE="${EVIDENCE}Even_coverage:${START_END_RATIO};"
  elif (( $(awk "BEGIN {print ($START_END_RATIO > 0.7 && $START_END_RATIO < 1.3)}") )); then
    CIRC_SCORE=$((CIRC_SCORE + 20))
    EVIDENCE="${EVIDENCE}Fair_coverage:${START_END_RATIO};"
  fi
  
  # Add score for split reads (up to 20 points)
  if [ "$SPLIT_COUNT" -ge 5 ]; then
    CIRC_SCORE=$((CIRC_SCORE + 20))
    EVIDENCE="${EVIDENCE}Split_reads:${SPLIT_COUNT}"
  elif [ "$SPLIT_COUNT" -ge 1 ]; then
    CIRC_SCORE=$((CIRC_SCORE + 10))
    EVIDENCE="${EVIDENCE}Split_reads:${SPLIT_COUNT}"
  fi
  
    # Determine status
    if [ "$CIRC_SCORE" -ge 60 ]; then  # Changed from 70
    STATUS="circular"
    elif [ "$CIRC_SCORE" -ge 30 ]; then  # Changed from 40
    STATUS="possibly_circular"
    else
    STATUS="linear"
    fi
  
  # Output result
  echo -e "$ref\t$CIRC_SCORE\t$STATUS\t$EVIDENCE" >> "$OUTPUT_TSV"
  
  # Clean up temp files
  rm "$START_READS" "$END_READS" "$SPANNING_READS"
  
done < "$SEQLEN_TEMP"

echo "Circularity analysis completed"

echo "DEBUG: $ref START_COV=$START_COV END_COV=$END_COV MIDDLE_COV=$MIDDLE_COV RATIO=$START_END_RATIO"

# Cleanup
rm "$COVERAGE_TEMP" "$SEQLEN_TEMP" "$DEPTH_ALL_TEMP"

echo "Mapping statistics written to $OUTPUT_TSV"