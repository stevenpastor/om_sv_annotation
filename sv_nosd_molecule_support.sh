#!/usr/bin/env bash

# download the segmental duplications track BED file for hg38 from UCSC:
echo "Downloading hg38 segmental duplications track from UCSC..."
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz && zcat genomicSuperDups.txt.gz | cut -f2-4 > genomicSuperDups.txt && rm genomicSuperDups.txt.gz

# sort and merge the segmental duplications BED file if within 50kbp:
echo ""
echo ""
echo "Sorting and merging hg38 segmental duplications track..."
#sort -k1,1 -k2,2n genomicSuperDups.txt | bedtools merge -d 50000 -i - > genomicSuperDups.sorted.merged.txt
sort -k1,1 -k2,2n genomicSuperDups.txt | bedtools merge -d 1 -i - > genomicSuperDups.sorted.merged.txt

# filter for >2kbp SVs and cut necessary cols:
echo ""
echo ""
echo "Filtering OM SVs >2kbp and replacing 23 and 24 for chrX and Y..."
grep -v "#" exp_refineFinal1_merged_filter.smap | awk '(($8-$7)-($6-$5))>2000 || (($8-$7)-($6-$5))<-2000' | awk '$8 != "-1.0"' | grep -v "trans" | cut -f1,2,3,5,6,7,8,10,18 | awk -v OFS='\t' '{ sub(/23$/, "X", $3) }1' | awk -v OFS='\t' '{ sub(/24$/, "Y", $3) }1' | awk '{ print $3"\t"$6"\t"$7"\t"$1","$2","$4","$5","$8","$9}' | awk -v OFS='\t' '{ sub(/\..$/, "", $2) }1' | awk -v OFS='\t' '{ sub(/\..$/, "", $3) }1' | sed 's/^/chr/g' | sort -k1,1 -k2,2n > optical_mapping_filtered_svs.txt

# obtain SVs not in SDs:
echo ""
echo ""
echo "Obtaining filtered OM SVs not in segmental duplications..."
bedtools intersect -a optical_mapping_filtered_svs.txt -b genomicSuperDups.sorted.merged.txt -v > optical_mapping_filtered_svs_not_in_sds.txt

echo ""
echo ""
echo "Initializing operations for OM SVs not in segmental duplications..."
# Counts number of molecules crossing an SV:
echo ""
echo ""
echo "Counting frequencies of molecules crossing each SV..."
while read line1 line2 line3 line4; do id=$(echo $line4 | cut -d',' -f2); start=$(echo $line4 | cut -d',' -f3); end=$(echo $line4 | cut -d',' -f4); awk -v i="$id" -v s="$start" -v e="$end" '{ if(s >= $6 && s <= $7) {print $0} else if(e >= $6 && e <= $7) {print $0} }' ../../../exp_refineFinal1_contig"$id".xmap | wc -l | sed 's/\ //g'; done <optical_mapping_filtered_svs_not_in_sds.txt > optical_mapping_filtered_svs_not_in_sds_number_molecules.txt

# Confidence of molecules crossing into the SV:
echo ""
echo ""
echo "Obtaining confidences of molecules crossing each SV..."
while read line1 line2 line3 line4; do id=$(echo $line4 | cut -d',' -f2); start=$(echo $line4 | cut -d',' -f3); end=$(echo $line4 | cut -d',' -f4); awk -v i="$id" -v s="$start" -v e="$end" '{ if(s >= $6 && s <= $7) {print $9} else if(e >= $6 && e <= $7) {print $9} }' ../../../exp_refineFinal1_contig"$id".xmap | sed 's/\ //g' | tr '\n' ',' | sed 's/,$//g'; echo ""; done <optical_mapping_filtered_svs_not_in_sds.txt > optical_mapping_filtered_svs_not_in_sds_confidence_only.txt

# Counts 5’ labels flanking an SV:
#In this file: optical_mapping_filtered_svs_not_in_sds.txt, want the 4th column but need the 2nd element (comma-delimited) in this list = contig ID:
echo ""
echo ""
echo "Counting frequency of labels before/after each molecule label directly flanking each SV..."
cut -f4 optical_mapping_filtered_svs_not_in_sds.txt | cut -d',' -f2,3,4 | tr ',' '\t' > tmp_contig_ids
#Use that contig ID to find it in the 1st column in this file: exp_refineFinal1_merged_q.cmap
#Once get the match (multiple rows), need to awk the 6th column in the QCMAP for the 3rd element in list (one output) and the 4th element (another output) and get only the 4th column in each instance:
while read line1 line2 line3; do awk -v i="$line1" -v s="$line2" -v e="$line3" 'function abs(x){return ((x < 0.0) ? -x : x)} {if ($1==i && abs($6-s) < 1) {print i","s","e","$4} else if($1==i && abs($6-e) < 1) {print i","s","e","$4}}' exp_refineFinal1_merged_q.cmap | tr '\n' '\t'; echo ""; done <tmp_contig_ids > tmp_contig_ids_labels

# Before labels 1 and 2 and after labels 1 and 2:
while read line1 line2; do id=$(echo $line1 | cut -d',' -f1); label1=$(echo $line1 | cut -d',' -f4); label2=$(echo $line2 | cut -d',' -f4); grep -v "#" ../../../exp_refineFinal1_contig"$id".xmap | cut -f14 | grep "($label1," | while read line; do echo $line | sed 's/(//g' | tr ')' '\n' | awk 'NF > 0' | cut -d',' -f1 | sed -n "0,/$label1/p" | wc -l; done < "${1:-/dev/stdin}" | awk '{$1=$1};1' | tr '\n' ',' | sed 's/,$//g'; echo "" | tr '\n' '\t'; grep -v "#" ../../../exp_refineFinal1_contig"$id".xmap | cut -f14 | grep "($label2," | while read line; do echo $line | sed 's/(//g' | tr ')' '\n' | awk 'NF > 0' | cut -d',' -f1 | sed -n "0,/$label2/p" | wc -l; done < "${1:-/dev/stdin}" | awk '{$1=$1};1' | tr '\n' ',' | sed 's/,$//g'; echo "" | tr '\n' '\t'; grep -v "#" ../../../exp_refineFinal1_contig"$id".xmap | cut -f14 | grep "($label1," | while read line; do echo $line | sed 's/(//g' | tr ')' '\n' | awk 'NF > 0' | cut -d',' -f1 | sed "1,/$label1/d" | wc -l; done < "${1:-/dev/stdin}" | awk '{$1=$1};1' | tr '\n' ',' | sed 's/,$//g'; echo "" | tr '\n' '\t'; grep -v "#" ../../../exp_refineFinal1_contig"$id".xmap | cut -f14 | grep "($label2," | while read line; do echo $line | sed 's/(//g' | tr ')' '\n' | awk 'NF > 0' | cut -d',' -f1 | sed "1,/$label2/d" | wc -l; done < "${1:-/dev/stdin}" | awk '{$1=$1};1' | tr '\n' ',' | sed 's/,$//g'; echo ""; done <tmp_contig_ids_labels > tmp_contig_ids_labels_before_after_label1_and_label2

# Paste together all files:
paste optical_mapping_filtered_svs_not_in_sds.txt optical_mapping_filtered_svs_not_in_sds_number_molecules.txt optical_mapping_filtered_svs_not_in_sds_confidence_only.txt tmp_contig_ids_labels_before_after_label1_and_label2 > svs_not_in_sds_support.txt

# Separate into 4 categories of molecule support confidence:
#1. molecule support >30, 5prime labels >=10, 3prime labels >=10, and average confidence >30:
while read line1 line2 line3 line4 line5 line6 line7 line8 line9 line10; do average=$(echo $line6 | tr ',' '\n' | awk '{ total += $1; count++ } END { print total/count }'); five=$(echo $line7 | tr ',' '\n' | awk '$1>=10' | wc -l); three=$(echo $line10 | tr ',' '\n' | awk '$1>=10' | wc -l); echo $line5";"$five";"$three";"$average";"$line1";"$line2";"$line3";"$line4";"$line6";"$line7";"$line8";"$line9";"$line10 | tr ';' '\t'; done <svs_not_in_sds_support.txt | awk -v a="$line1" -v b="$line2" -v c="$line3" -v d="$line4" -v e="$line5" -v f="$line6" -v g="$line7" -v h="$line8" -v i="$line9" -v j="$line10" '$1>30 && $2>=10 && $3>=10 && $4>30 {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",a,b,c,d,e,f,g,h,i,j}' > svs_not_in_sds_support_highestconf.txt
#2. molecule support >10 but <30:
while read line1 line2 line3 line4 line5 line6 line7 line8 line9 line10; do average=$(echo $line6 | tr ',' '\n' | awk '{ total += $1; count++ } END { print total/count }'); five=$(echo $line7 | tr ',' '\n' | awk '$1>=10' | wc -l); three=$(echo $line10 | tr ',' '\n' | awk '$1>=10' | wc -l); echo $line5";"$five";"$three";"$average";"$line1";"$line2";"$line3";"$line4";"$line6";"$line7";"$line8";"$line9";"$line10 | tr ';' '\t'; done <svs_not_in_sds_support.txt | awk -v a="$line1" -v b="$line2" -v c="$line3" -v d="$line4" -v e="$line5" -v f="$line6" -v g="$line7" -v h="$line8" -v i="$line9" -v j="$line10" '$1>10 && $1<30 {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",a,b,c,d,e,f,g,h,i,j}' > svs_not_in_sds_support_moderateconf.txt
# 3. lowest molecule support:
while read line1 line2 line3 line4 line5 line6 line7 line8 line9 line10; do average=$(echo $line6 | tr ',' '\n' | awk '{ total += $1; count++ } END { print total/count }'); five=$(echo $line7 | tr ',' '\n' | awk '$1>=10' | wc -l); three=$(echo $line10 | tr ',' '\n' | awk '$1>=10' | wc -l); echo $line5";"$five";"$three";"$average";"$line1";"$line2";"$line3";"$line4";"$line6";"$line7";"$line8";"$line9";"$line10 | tr ';' '\t'; done <svs_not_in_sds_support.txt | awk -v a="$line1" -v b="$line2" -v c="$line3" -v d="$line4" -v e="$line5" -v f="$line6" -v g="$line7" -v h="$line8" -v i="$line9" -v j="$line10" '$1<10 {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",a,b,c,d,e,f,g,h,i,j}' > svs_not_in_sds_support_lowestconf.txt

# Clean-up:
rm -f optical_mapping_filtered_svs_not_in_sds_number_molecules.txt; rm -f optical_mapping_filtered_svs_not_in_sds_confidence_only.txt; rm -f tmp_contig_ids_labels_before_after_label1_and_label2

