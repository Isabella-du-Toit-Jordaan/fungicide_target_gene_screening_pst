#!/bin/bash
set -euo pipefail
shopt -s nullglob

WORKDIR="/home/2016042547/FUS_Sequencing/2023+2025_FUS_exonic_variants"
cd "$WORKDIR"

READS_DIR="reads"
REF="ref/genes_combined.fasta"
EXON_TAB="exons.tab"   # CHROM  START(1-based) END(1-based)
THREADS=8
GENES=(cyp51 cytb-1 cytb-2 cytb-3 sdhb sdhc sdhd)

# output dirs
mkdir -p align medaka_var medaka_vcf mpileup_vcf renamed_vcf annotated_vcf per_gene tmp_sites \
         qc_stats logs tmp_work filtered_vcf genotype_tables genotypes_variants genotypes_nonvariant genotypes_all

# load modules
module load minimap2
module load samtools
module load medaka
module load bcftools

echo "[INFO] Starting pipeline in $PWD"

if [ ! -f "${REF}.fai" ]; then
    echo "[INFO] Creating FASTA index for $REF"
    samtools faidx "$REF"
fi

# Step 1: index reference

echo "[STEP 1] Indexing reference"
minimap2 -d ref/genes_combined.mmi "$REF"


# Step 2: mapping -> BAM

echo "[STEP 2] Mapping reads to BAM"
for fq in "${READS_DIR}"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    echo "[MAP] $sample"
    minimap2 -t "$THREADS" -ax map-ont ref/genes_combined.mmi "$fq" \
        | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "align/${sample}.sorted.bam" -
    samtools index "align/${sample}.sorted.bam"
done


# Step 3: quick QC

echo "[STEP 3] Quick QC"
mkdir -p qc_stats
for bam in align/*.sorted.bam; do
    sample=$(basename "$bam" .sorted.bam)
    samtools flagstat "$bam" > qc_stats/"$sample".flagstat
    samtools idxstats "$bam" > qc_stats/"$sample".idxstats
done


# Step 4: split BAMs by gene

echo "[STEP 4] Splitting BAMs by gene"
for bam in align/*.sorted.bam; do
    sample=$(basename "$bam" .sorted.bam)
    for g in "${GENES[@]}"; do
        outbam="align/${sample}_${g}.bam"
        echo "[SPLIT] $sample -> $g"
        samtools view -b "$bam" "$g" > "$outbam"
        samtools index "$outbam"
    done
done


# Step 4.5: Calculate coverage stats

echo "[STEP 4.5] Calculating coverage stats"
mkdir -p coverage_stats
echo -e "Sample\tGene\tReadsMapped\tMedianReadLength" > coverage_stats/coverage_summary.tsv
for bam in align/*_*.bam; do
    base=$(basename "$bam" .bam)
    sample=${base%%_*}
    gene=${base#*_}
    mapped=$(samtools view -c -F 260 "$bam")
    if [ "$mapped" -gt 0 ]; then
        median_len=$(samtools view "$bam" | awk '{print length($10)}' | sort -n | \
            awk '{a[NR]=$1} END {if (NR==0) print 0; else {mid=int((NR+1)/2); if (NR%2) print a[mid]; else print (a[mid]+a[mid+1])/2}}')
    else
        median_len=0
    fi
    echo -e "${sample}\t${gene}\t${mapped}\t${median_len}" >> coverage_stats/coverage_summary.tsv
done


# Step 5: Medaka variant calling

echo "[STEP 5] Running Medaka variant calling"
for bam in align/*_*.bam; do
    base=$(basename "$bam" .bam)   # sample01_cyp51
    outdir="medaka_var/${base}"
    mkdir -p "$outdir"
    echo "[MEDAKA] $base"
    medaka_variant -i "$bam" -r "$REF" -o "$outdir" -t "$THREADS"
    if [ ! -f "$outdir/medaka.vcf" ]; then
        echo "[WARN] Missing medaka.vcf for $base; skipping"
        continue
    fi
    bgzip -c "$outdir/medaka.vcf" > "medaka_vcf/${base}.vcf.gz"
    bcftools index -c "medaka_vcf/${base}.vcf.gz"
done


# Step 6: Medaka sites to per-base BED

echo "[STEP 6] Extracting Medaka-called positions to BED"
mkdir -p tmp_sites
for v in medaka_vcf/*.vcf.gz; do
    base=$(basename "$v" .vcf.gz)
    bcftools query -f '%CHROM\t%POS0\t%POS\n' "$v" > "tmp_sites/${base}.bed" || true
done


# Step 7: Run mpileup at EXON positions (using exons.tab)

echo "[STEP 7] Running mpileup to call variants inside exon regions"
mkdir -p mpileup_vcf tmp_work

# Create per-gene 0-based exon BEDs in tmp_work
rm -f tmp_work/*.exons.bed
while IFS=$'\t' read -r chrom start end; do
    start0=$((start-1))
    echo -e "${chrom}\t${start0}\t${end}" >> tmp_work/${chrom}.exons.bed
done < "$EXON_TAB"

# For every sample_gene BAM, run mpileup restricted to gene exons
for bam in align/*_*.bam; do
    base=$(basename "$bam" .bam)    # sample01_cyp51
    gene=${base#*_}                 # cyp51
    exon_bed="tmp_work/${gene}.exons.bed"
    out="mpileup_vcf/${base}.vcf.gz"
    echo "[MPILEUP] $base (exons: $exon_bed)"
    if [ ! -s "$exon_bed" ]; then
        echo "[MPILEUP] No exon intervals for $gene -> skipping mpileup for $base"
        continue
    fi
    bcftools mpileup -Ou -f "$REF" -R "$exon_bed" "$bam" -a FORMAT/AD,FORMAT/DP \
      | bcftools call -mv -Oz -o "$out"
    bcftools index -c -f "$out" || true
done


# Step 8: Reheader and annotate Medaka VCF with AD/DP from mpileup

echo "[STEP 8] Reheadering & annotating Medaka VCFs with AD/DP from mpileup"
mkdir -p renamed_vcf annotated_vcf

for medaka in medaka_vcf/*.vcf.gz; do
    base=$(basename "$medaka" .vcf.gz)   # sample01_cyp51
    mpileup_v="mpileup_vcf/${base}.vcf.gz"
    medaka_re="renamed_vcf/${base}.medaka.rehead.vcf.gz"
    mpileup_re="tmp_work/${base}.mpileup.rehead.vcf.gz"
    annotated_out="annotated_vcf/${base}.vcf.gz"

    echo "[ANNOTATE] $base"

    if [ ! -f "$mpileup_v" ]; then
        echo "[WARN] No mpileup VCF ($mpileup_v) for $base -> create empty (no records) placeholder"
        bcftools view -h "$medaka" | bgzip -c > "$mpileup_v" || true
        touch "${mpileup_v}.csi" || true
    fi

    echo "[REHEADER] mpileup -> $mpileup_re (sample=$base)"
    bcftools reheader -s <(echo "$base") -o "$mpileup_re" "$mpileup_v" || true
    bcftools index -c -f "$mpileup_re" || true

    echo "[REHEADER] medaka -> $medaka_re (sample=$base)"
    bcftools reheader -s <(echo "$base") -o "$medaka_re" "$medaka" || true
    bcftools index -c -f "$medaka_re" || true

    medaka_sample=$(bcftools query -l "$medaka_re" 2>/dev/null || echo "")
    mpileup_sample=$(bcftools query -l "$mpileup_re" 2>/dev/null || echo "")
    if [[ -z "$medaka_sample" || -z "$mpileup_sample" || "$medaka_sample" != "$mpileup_sample" ]]; then
        echo "[ERROR] sample name mismatch after reheader: medaka='$medaka_sample' mpileup='$mpileup_sample' -> skipping"
        rm -f "$medaka_re" "$mpileup_re"
        continue
    fi

    echo "[ANNOTATE] Adding FORMAT/AD,FORMAT/DP to $medaka_re -> $annotated_out"
    bcftools annotate -a "$mpileup_re" -c CHROM,POS,REF,ALT,FORMAT/AD,FORMAT/DP -Oz -o "$annotated_out" "$medaka_re" || true
    bcftools index -c -f "$annotated_out" || true

    rm -f "$medaka_re" "$mpileup_re"
done


# Step 8.5: Filter annotated VCFs

echo "[STEP 8.5] Filtering annotated VCFs (GQ>=20, DP>=50, AD>=15, AF thresholds)"
mkdir -p filtered_vcf

for vcf in annotated_vcf/*.vcf.gz; do
    base=$(basename "$vcf" .vcf.gz)
    echo "[FILTER] $base"

    tmp_tags="tmp_work/${base}.tags.vcf.gz"
    bcftools +fill-tags "$vcf" -Oz -o "$tmp_tags" -- -t AF || true
    bcftools index -c -f "$tmp_tags" || true

    out="filtered_vcf/${base}.filtered_snps.vcf.gz"
    # SNP filter
    bcftools view -i 'QUAL>=5 && FORMAT/GQ>=20 && FORMAT/DP>=50 && FORMAT/AD[*:1]>=15 && ((INFO/AF>=0.30 && INFO/AF<=0.70) || INFO/AF>=0.85)' \
        -Oz -o "$out" "$tmp_tags" || true

    # INDELs (kept as separate file)
    bcftools view -i 'TYPE="indel" && FORMAT/DP>=100 && FORMAT/AD[*:1]>=40 && INFO/AF>=0.90' \
        -Oz -o "${out%.vcf.gz}.indels.vcf.gz" "$tmp_tags" || true

    bcftools index -c -f "$out" || true

    rm -f "$tmp_tags" "$tmp_tags.csi" || true
done


# Step 9: Build per-gene matrices

echo "[STEP 9] Building per-gene matrices"
mkdir -p per_gene
for gene in "${GENES[@]}"; do
    echo "[MATRIX] gene=$gene"
    out="per_gene/${gene}.tsv"
    echo -e "Sample\tCHROM\tPOS\tREF\tALT\tGT\tAD\tDP\tAF" > "$out"
    awk -v g="$gene" '$1==g' "$EXON_TAB" > "tmp_work/${gene}.bed"
    for v in filtered_vcf/*.vcf.gz; do
        base=$(basename "$v" .vcf.gz)
        bcftools view -R "tmp_work/${gene}.bed" "$v" \
          | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD\t%DP]\n' \
          | awk -v S="$base" 'BEGIN{OFS="\t"}{
              gt=$5; ad=$6; dp=$7;
              if(ad=="." || dp==".") next;
              split(ad,a,","); ref=a[1]+0; alt=a[2]+0;
              af = alt/(ref+alt+0.001);
              print S,$1,$2,$3,$4,gt,ad,dp,af
            }' >> "$out" || true
    done
done


# Step 10: Summarize per-gene matrices for downstream analysis (tidy tables for R analysis)

echo "[STEP10] Starting tidy table generation"
mkdir -p tidy_tables genotypes_variants genotypes_nonvariant genotypes_all tmp_work

# ensure fasta index exists (used to fetch reference base)
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

for gene in "${GENES[@]}"; do
    echo "[STEP10] Processing ${gene}"

    # find sample VCFs for this gene (ignore .indels)
    mapfile -t SAMPLE_VCFS < <(ls -1 filtered_vcf/sample*_"${gene}".filtered_snps.vcf.gz 2>/dev/null || true)

    # Build list of sample names (strip suffix _${gene}.filtered_snps.vcf.gz)
    SAMPLE_NAMES=()
    for f in "${SAMPLE_VCFS[@]}"; do
        name=$(basename "$f")
        samp=$(printf "%s" "$name" | sed -E "s/_${gene}\.filtered_snps\.vcf\.gz$//")
        SAMPLE_NAMES+=("$samp")
    done

    # Variants
    VAR_OUT="tidy_tables/${gene}.genotypes_variants.tsv"
    echo -e "sample\tgene\tchrom\tpos\tref\talt\tgt_score\tGT\tAD\tDP" > "$VAR_OUT"

    if [ ${#SAMPLE_VCFS[@]} -eq 0 ]; then
        echo "  No sample VCFs found for ${gene}"
        : > "$VAR_OUT"
    else
        for vcf in "${SAMPLE_VCFS[@]}"; do
            sample=$(basename "$vcf" | sed -E "s/_${gene}\.filtered_snps\.vcf\.gz$//")

            if [ "$(bcftools view -H "$vcf" 2>/dev/null | wc -l)" -eq 0 ]; then
                continue
            fi

            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%DP]\n' "$vcf" 2>/dev/null \
            | awk -v SAMP="$sample" -v GENE="$gene" 'BEGIN{OFS="\t"}
                {
                  chrom=$1; pos=$2; ref=$3; alt=$4; gt=$5; ad=$6; dp=$7;

                  # missing
                  if(gt==".") {
                      print SAMP,GENE,chrom,pos,ref,alt,"NA",gt,ad,dp;
                      next;
                  }

                  # haploid GT (Medaka 0 or 1) -> infer using AD
                  if(gt ~ /^[0-9]+$/) {
                      split(ad,a,",");
                      refc=a[1]+0; altc=a[2]+0;
                      total=refc+altc+0.001;
                      af=altc/total;

                      if(af < 0.2) code=0;
                      else if(af >= 0.3 && af <= 0.7) code=0.5;
                      else if(af >= 0.85) code=1;
                      else code="NA";

                      print SAMP,GENE,chrom,pos,ref,alt,code,gt,ad,dp;
                      next;
                  }

                  # diploid/phased GT (0/0, 0/1, 1/1)
                  n=split(gt,a,/[/|]/); miss=0; altcount=0;
                  for(i=1;i<=n;i++){
                    if(a[i]==".") miss=1;
                    else if(a[i]!="0") altcount++;
                  }

                  if(miss) code="NA";
                  else {
                    frac = altcount / n;
                    if(frac < 0.2) code=0;
                    else if(frac >= 0.3 && frac <= 0.7) code=0.5;
                    else if(frac >= 0.85) code=1;
                    else code="NA";
                  }

                  print SAMP,GENE,chrom,pos,ref,alt,code,gt,ad,dp;
                }' >> "$VAR_OUT"
        done
    fi

    # all positions
    ALL_POS_FILE="tmp_work/${gene}.all_pos.txt"
    awk -v g="$gene" '$1==g { for(i=$2;i<=$3;i++) print $1"\t"i }' "$EXON_TAB" \
        | sort -k1,1 -k2,2n > "$ALL_POS_FILE"

    if [ ! -s "$ALL_POS_FILE" ]; then
        echo "  No exons found for ${gene}"
        : > "$ALL_POS_FILE"
    fi

    ALL_OUT="tidy_tables/${gene}.genotypes_all.tsv"
    NONVAR_OUT="tidy_tables/${gene}.genotypes_nonvariants.tsv"

    echo -e "sample\tgene\tchrom\tpos\tref\talt\tgt_score\tGT\tAD\tDP" > "$ALL_OUT"
    echo -e "sample\tgene\tchrom\tpos\tref\talt\tgt_score\tGT\tAD\tDP" > "$NONVAR_OUT"

    BASE_FILE="tmp_work/${gene}.all_base.tsv"
    : > "$BASE_FILE"

    while IFS=$'\t' read -r chrom pos; do
        refbase=$(samtools faidx "$REF" "${chrom}:${pos}-${pos}" 2>/dev/null | sed -n '2p' || echo ".")
        for samp in "${SAMPLE_NAMES[@]}"; do
            printf "%s\t%s\t%s\t%s\t%s\t.\t0\t.\t.\t.\n" \
                "$samp" "$gene" "$chrom" "$pos" "$refbase" >> "$BASE_FILE"
        done
    done < "$ALL_POS_FILE"

    if [ ${#SAMPLE_NAMES[@]} -eq 0 ]; then
        echo "  No samples for ${gene}"
        continue
    fi

    VAR_LOOKUP="tmp_work/${gene}.var_lookup.tsv"
    awk -v OFS="\t" 'NR>1 {
        key=$1"\t"$2"\t"$3"\t"$4;
        payload=$5"\t"$6"\t"$7"\t"$8"\t"$9;
        print key,payload
    }' "$VAR_OUT" > "$VAR_LOOKUP" || true

    awk -v VARL="$VAR_LOOKUP" 'BEGIN{
         while((getline < VARL) > 0){
           key=$1"\t"$2"\t"$3"\t"$4;
           payload=$5"\t"$6"\t"$7"\t"$8"\t"$9;
           V[key]=payload;
         }
         close(VARL);
       }
       {
         key=$1"\t"$2"\t"$3"\t"$4;
         if(key in V){
            split(V[key],P,"\t");
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               $1,$2,$3,$4,P[1],P[2],P[3],P[4],P[5],P[6];
         } else {
            printf "%s\t%s\t%s\t%s\t%s\t.\t0\t.\t.\t.\n",
               $1,$2,$3,$4,$5;
         }
       }' "$BASE_FILE" >> "$ALL_OUT"

    awk 'NR>1 && $7==0 { print }' "$ALL_OUT" >> "$NONVAR_OUT"

    # legacy variants
    LEG_VARIANTS="genotypes_variants/${gene}_genotypes_variants.tsv"
    : > "$LEG_VARIANTS"

    awk 'NR>1 { print $3"\t"$4"\t"$5"\t"$6 }' "$VAR_OUT" | sort -u \
    | while IFS=$'\t' read -r chrom pos ref alt; do
        line="${chrom}\t${pos}\t${ref}\t${alt}"
        for samp in "${SAMPLE_NAMES[@]}"; do
            gt=$(awk -v s="$samp" -v c="$chrom" -v p="$pos" \
                  'BEGIN{FS="\t"} $1==s && $3==c && $4==p {print $8; exit}' "$VAR_OUT")
            if [ -z "$gt" ]; then gt="."; fi
            line="${line}\t${samp}=${gt}"
        done
        echo -e "$line" >> "$LEG_VARIANTS"
    done

    cp "$LEG_VARIANTS" "genotypes_all/${gene}_genotypes_all.tsv" || true
    cp "$NONVAR_OUT" "genotypes_nonvariant/${gene}_genotypes_nonvariant.tsv" || true

    rm -f "$VAR_LOOKUP" "$BASE_FILE" "$ALL_POS_FILE"
done

echo "[STEP10] Finished tidy table generation"
