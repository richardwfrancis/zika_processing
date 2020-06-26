## zika_processing
important note that vt has been altered in ~/software/vt
see https://github.com/atks/vt/issues/46 for the change required

original vcf data obtained from Alexia. This had been run through GATK and samtools and the product is the intersect from those tools. This data is in ./raw
this vcf file and an ANNOVAR annotated version are in the folder raw_all_vcf
there are 660,007 variants in these files.
```{bash}
ls raw
```

variant data had not been decomposed or normalised so this was performed via the seqnextgen pipeline, see sng_script.sh
```{bash}
cat sng_script.sh
```

this copies the original vcf into ./tmp and normalises and decomposes with vt
the normalised VCF is then annotated with VEP
the final file combined_HaplotypeCaller.d.n.vep.vcf.gz is a bgzipped, tabix and grabix indexed VCF containing 661866 variants (1859 more than the original due to the normalisation process)
note that the data is loaded to GEMINI with the --passonly flag and so is equivalent to the combined_HaplotypeCaller.d.n.vep.PASSonly.vcf file compared to dbNSFP below
```{bash}
ls tmp
```

filter the VCF to PASS only
this file conains variants that were deemed false positives by GATK (see VQSR docs https://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr)
it was agreed that this data should only contain PASSed variants so filter this data accordingly and index
combined_HaplotypeCaller.d.n.vep.PASSonly.vcf contains 612761 variants
```{bash}
bcftools view -f PASS combined_HaplotypeCaller.d.n.vep.vcf > combined_HaplotypeCaller.d.n.vep.PASSonly.vcf
bgzip combined_HaplotypeCaller.d.n.vep.PASSonly.vcf 
tabix combined_HaplotypeCaller.d.n.vep.PASSonly.vcf.gz
grabix index combined_HaplotypeCaller.d.n.vep.PASSonly.vcf.gz
```

compare this file to dbNSFP to find ns and sv variants
note that it is not clear whether or not dbNSFP is normalised and decomposed however doing this with a non-normalised file matched 110 fewer variants
the dbNSFP data is installed in /data/workspace/alexia/dbNSFP and the following was executed in this directory to allow access ot the dbNSFP suite 
```{bash}
java -Xmx5g search_dbNSFP40a -i /data/workspace/richard/zika/tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf -o /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out -v hg19 > /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out.log 2>&1
```

the output of the above produced
```
58111 SNP(s) are found. Written to /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out
554650 SNP(s) are not found. Written to /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out.err
```
however the above PASSonly.out file contains 58556 lines of variant data (445 more than stated due to duplicates (e.g. 19:44919189 and 16:1992649) arising from alternative consequences in different transcripts)

to get the required annotation for input into the acmg tool the VCF was brought into GEMINI
relevant data fields were manually extracted from within the resulting sqlite database
```{sql}
sqlite3 Illumina.db
.headers on
.mode tabs
.output zika_gemini_export.tsv
select chrom,start,end,ref,alt,clinvar_sig,impact,gene,exon,vep_hgvsc,vep_hgvsp,aa_change,codon_change,pfam_domain,max_aaf_all,polyphen_score,sift_score,cadd_scaled from variants;
```

find exported GEMINI data that was found in dbNSFP
use GEMINI chr start + 1 == sbNSFP hg19_pos(1-based)
dbnsfp has 470 columns so just narrow it down to the hg19 columns we want. include the ref and alt to make sure we are matching the right variant found in dbNSFP.
"ref" "alt" "hg19_chr" "hg19_pos(1-based)"
```{bash}
awk 'BEGIN{FS="\t"; OFS="\t";}{print $3,$4,$8,$9; }' zika_vep_dbnsfp.PASSonly.out > zika_vep_dbnsfp.PASSonly.hg19cols.out
```

use R to do the mergin and output a file for input into acmg
this file contains 58556 lines of variant data representing the subset of the GEMINI data that also matched to dbNSFP.
note that this includes the 445 duplicates mentioned above
```{r}
library("tidyverse")
dbnsfp <- read_tsv("zika_vep_dbnsfp.PASSonly.hg19cols.out", col_types="cccn")
gemini <- read_tsv("zika_gemini_export.tsv", col_types = "ciicccccccccccnnnn")
gemini <- gemini %>% mutate(chrom=gsub("chr","",chrom), newstart=start+1)
combined4 <- dbnsfp %>%
 left_join(gemini, by=c("hg19_chr"="chrom","hg19_pos(1-based)"="newstart","ref"="ref","alt"="alt"), keep=TRUE) %>%
 select('hg19_chr','start','hg19_pos(1-based)','end','ref','alt','clinvar_sig','impact','gene','exon','vep_hgvsc','vep_hgvsp','aa_change','codon_change','pfam_domain','max_aaf_all','polyphen_score','sift_score','cadd_scaled') %>%
 rename('chrom'='hg19_chr','newstart'='hg19_pos(1-based)')

write_tsv(combined4,"combined_dbnsfp_geminiPASSonly.tsv")
```

run ACMG
there are 58556 lines in the input file which will take a while to run
as there is interaction with live databases across networks, things can go wrong
to save re-running all 60K variants the file was split into 30 files containing roughly 2K variants each and a job script made
once completed combine the files
these files are in ./splitfiles
note there were 6 variants that had connection errors during the process. Their data (see linesWithConnectionErrors.txt) were manually replaced in all_acmg.out.tsv folloing re-running with acmg
```{bash}
split -l 2000 combined_dbnsfp_geminiPASSonly.tsv acmg_job_
for i in acmg*; do sed -i '1s/^/chrom\tstart\tnewstart\tend\tref\talt\tclinvar_sig\timpact\tgene\texon\tvep_hgvsc\tvep_hgvsp\taa_change\tcodon_change\tpfam_domain\tmax_aaf_all\tpolyphen_score\tsift_score\tcadd_scaled\n/' $i; done
for i in acmg*; do echo "python3 ../../software/ppTimo/Phenoparser/scripts/acmg.py -s $i -p $i -g /data/scratch/reference/acmg/acmg_id.db > $i"_acmg.log" 2>&1"; done > runacmg.sh
chmod +x runacmg.sh
./runacmg.sh
head -n 1 acmg_job_aaacmg.out.tsv > combined_dbnsfp_geminiPASSonly_acmg.tsv
awk 'BEGIN{FS="\t";} FNR>1 {print;}' acmg_job_*acmg.out.tsv >> combined_dbnsfp_geminiPASSonly_acmg.tsv
```

make a new VCF file of this subset of variants to annotate
first use zika_vep_dbnsfp.PASSonly.hg19cols.out to make a file to search for these variants in the original VCF
then grep out the header lines and dbNSFP variants
note that this file contains 58111 variants which excludes the 445 duplicate findings in dbNSFP mentioned above
```{bash}
awk 'BEGIN{FS="\t";OFS="\t";} FNR>1 {print "chr"$3,$4,".",$1,$2; }' zika_vep_dbnsfp.PASSonly.hg19cols.out > getdbNSFPsubset.txt
grep -e "^#" tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf > ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf
grep -f getdbNSFPsubset.txt tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf >> ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf
```

 edit ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf to add the dbNSFP search command
```{bash}
vi ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf
##dbNSFPCommand="java -Xmx5g search_dbNSFP40a -i /data/workspace/richard/zika/tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf -o /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out -v hg19 > /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out.log 2>&1"
```

make a file containing the acmg annotation for these variants to add to the VCF annotation
```{bash}
awk 'BEGIN{FS="\t";OFS="\t";} FNR >1 {print "chr"$1,$3,$5,$6,$20"|"$22"|"$24"|"$26"|"$28"|"$30"|"$32"|"$34"|"$36; }' combined_dbnsfp_geminiPASSonly_acmg.tsv | sort -k1,1V -k2,2n > combined_dbnsfp_geminiPASSonly_acmg_anno.tsv
bgzip combined_dbnsfp_geminiPASSonly_acmg_anno.tsv
tabix -s1 -b2 -e2 combined_dbnsfp_geminiPASSonly_acmg_anno.tsv.gz
```

annotate dbNSFP variant VCF with ACMG
```{bash}
echo '##INFO=<ID=ACMG,Number=.,Type=String,Description="ACMG classification and evidence. Format: classification|pvs1|ps1|pm1|pm2|pm4|pm5|pp2|pp3">' > combined_dbnsfp_geminiPASSonly_acmg_anno.hdr.tsv
bcftools annotate -a combined_dbnsfp_geminiPASSonly_acmg_anno.tsv.gz -h combined_dbnsfp_geminiPASSonly_acmg_anno.hdr.tsv -c CHROM,POS,REF,ALT,ACMG -o ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf -O v ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf
bgzip ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf
tabix ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf.gz
```

#########################
4/6/20 - a couple of issues with the acmg code were spotted when responding to reviewers
a new set of LP and P classifications were made and need to be added to the VCF file
acmg_pvs1pm1Check_anno.txt - contains the correct LP and P data
combined_dbnsfp_geminiPASSonly_acmg.tsv - contains all the acmg results (also contains erroneous pm1 and pvs1 data)
I can correct the pm1 data by looking at the pfam_domain and domains columns
I cannot correct the pvs1 column without a fresh acmg database
So I have bitten the bullet and generated a new acmg database using the 242 genes from combined_dbnsfp_geminiPASSonly_acmg.tsv that were deemed pvs1
```{bash}
awk -v FS="\t" '{if ($22 == 1){ print $9; } }' combined_dbnsfp_geminiPASSonly_acmg.tsv | sort -u > combined_dbnsfp_geminiPASSonly_acmg_uniqueGenes_pvs1.tsv
# run on my computer in /Users/rfrancis_adm/Documents/TKI/GeneticsHealth/Timo/zika/scientific_data_paper/submitted_201219/responseToReviewers/files4table2/createNewacmgdb
python3 ../../Phenoparser/scripts/get_gene_disease.py -k <OMIM_KEY> -g combined_dbnsfp_geminiPASSonly_acmg_uniqueGenes_pvs1.tsv -t file -o acmg_id_zika1.db > gd1.log 2>&1
```
Need to use this to merge the disease data onto this file for each variant
Then need to reassign the 1s and 0s for pvs1 and pm1
	For any of the corrected data I need to reassign the acmg classifier 
All of the data in acmg_pvs1pm1Check_anno.txt is correct so use that in preference to the data in combined_dbnsfp_geminiPASSonly_acmg.tsv
Make a final corrected file to use in the annotation
Annotate as below
```{bash}
echo '##INFO=<ID=ACMG,Number=.,Type=String,Description="ACMG classification and evidence. Format: classification|pvs1|ps1|pm1|pm2|pm4|pm5|pp2|pp3">' > combined_dbnsfp_geminiPASSonly_acmg_anno.hdr.tsv
bcftools annotate -a combined_dbnsfp_geminiPASSonly_acmg_anno.tsv.gz -h combined_dbnsfp_geminiPASSonly_acmg_anno.hdr.tsv -c CHROM,POS,REF,ALT,ACMG -o ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf -O v ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf
bgzip ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf
tabix ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf.gz
```
#########################

filter 3 samples out of the final files
```{bash}
bcftools view -S ^samples2exclude.txt -O z -o combined_HaplotypeCaller.d.n.vep.PASSonly.45samples.vcf.gz combined_HaplotypeCaller.d.n.vep.PASSonly.vcf.gz
bcftools view -S ^samples2exclude.txt -O z -o combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg_45samples.vcf.gz combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf.gz
```

final vcf headers
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -4.7548 <= x < -0.2564">
##FILTER=<ID=VQSRTrancheINDEL99.90to100.00+,Description="Truth sensitivity tranche level for INDEL model at VQS Lod < -546.34">
##FILTER=<ID=VQSRTrancheINDEL99.90to100.00,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -546.34 <= x < -4.7548">
##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -8.7387 <= x < -0.7207">
##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description="Truth sensitivity tranche level for SNP model at VQS Lod < -8137.5268">
##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -8137.5268 <= x < -8.7387">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
##GATKCommandLine=<ID=ApplyVQSR,CommandLine="ApplyVQSR  --recal-file /SCRATCH/aweeks_analysis/ZIKA/variants/genotype_all.snp.recal.csv --tranches-file variants/genotype_all.snp.tranches --output variants/genotype_all.snp_recalibrated.vcf --truth-sensitivity-filter-level 99.0 --mode SNP --variant variants/genotype_all.vcf --intervals /SCRATCH/programs/misc/sureselect_v6UTR_bed/S07604624_Regions.bed --interval-padding 100 --reference /SCRATCH/reference/hg19.fa  --use-allele-specific-annotations false --ignore-all-filters false --exclude-filtered false --interval-set-rule UNION --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --disable-tool-default-read-filters false",Version=4.0.2.1,Date="March 29, 2019 7:07:01 AM AWST">
##GATKCommandLine=<ID=ApplyVQSR,CommandLine="ApplyVQSR  --recal-file /SCRATCH/aweeks_analysis/ZIKA/variants/genotype_all.snp_recalibrated.indel.recal.csv --tranches-file variants/genotype_all.snp_recalibrated.indel.tranches --output variants/genotype_all.snp_recalibrated.indel_recalibrated.vcf --truth-sensitivity-filter-level 99.0 --mode INDEL --variant variants/genotype_all.snp_recalibrated.vcf --intervals /SCRATCH/programs/misc/sureselect_v6UTR_bed/S07604624_Regions.bed --interval-padding 100 --reference /SCRATCH/reference/hg19.fa  --use-allele-specific-annotations false --ignore-all-filters false --exclude-filtered false --interval-set-rule UNION --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --disable-tool-default-read-filters false",Version=4.0.2.1,Date="March 29, 2019 7:10:58 AM AWST">
##GATKCommandLine=<ID=CombineGVCFs,CommandLine="CombineGVCFs  --output variants/combined.g.vcf --variant variants/variant.list --intervals /SCRATCH/programs/misc/sureselect_v6UTR_bed/S07604624_Regions.bed --interval-padding 100 --reference /SCRATCH/reference/hg19.fa  --annotation-group StandardAnnotation --disable-tool-default-annotations false --convert-to-base-pair-resolution false --break-bands-at-multiples-of 0 --ignore-variants-starting-outside-interval false --interval-set-rule UNION --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --disable-tool-default-read-filters false",Version=4.0.2.1,Date="March 29, 2019 3:04:22 AM AWST">
##GATKCommandLine=<ID=GenotypeGVCFs,CommandLine="GenotypeGVCFs  --output variants/genotype_all.vcf --max-alternate-alleles 40 --variant variants/combined.g.vcf --intervals /SCRATCH/programs/misc/sureselect_v6UTR_bed/S07604624_Regions.bed --interval-padding 100 --reference /SCRATCH/reference/hg19.fa  --use-new-qual-calculator false --annotate-with-num-discovered-alleles false --heterozygosity 0.001 --indel-heterozygosity 1.25E-4 --heterozygosity-stdev 0.01 --standard-min-confidence-threshold-for-calling 10.0 --max-genotype-count 1024 --sample-ploidy 2 --annotation-group StandardAnnotation --disable-tool-default-annotations false --only-output-calls-starting-in-intervals false --interval-set-rule UNION --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --disable-tool-default-read-filters false",Version=4.0.2.1,Date="March 29, 2019 5:39:29 AM AWST">
##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller  --dbsnp /SCRATCH/programs/gatk-4.0.2.1/hg19/dbsnp_138.hg19.vcf --emit-ref-confidence GVCF --max-alternate-alleles 50 --output variants/JN-163.dedup.raw.snps.indels.g.vcf --intervals /SCRATCH/programs/misc/sureselect_v6UTR_bed/S07604624_Regions.bed --interval-padding 100 --input alignments/JN-163.dedup.bam --reference /SCRATCH/reference/hg19.fa  --annotation-group StandardAnnotation --annotation-group StandardHCAnnotation --disable-tool-default-annotations false --gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 --gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 --gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 --gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 --gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30 --gvcf-gq-bands 31 --gvcf-gq-bands 32 --gvcf-gq-bands 33 --gvcf-gq-bands 34 --gvcf-gq-bands 35 --gvcf-gq-bands 36 --gvcf-gq-bands 37 --gvcf-gq-bands 38 --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 41 --gvcf-gq-bands 42 --gvcf-gq-bands 43 --gvcf-gq-bands 44 --gvcf-gq-bands 45 --gvcf-gq-bands 46 --gvcf-gq-bands 47 --gvcf-gq-bands 48 --gvcf-gq-bands 49 --gvcf-gq-bands 50 --gvcf-gq-bands 51 --gvcf-gq-bands 52 --gvcf-gq-bands 53 --gvcf-gq-bands 54 --gvcf-gq-bands 55 --gvcf-gq-bands 56 --gvcf-gq-bands 57 --gvcf-gq-bands 58 --gvcf-gq-bands 59 --gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99 --indel-size-to-eliminate-in-ref-model 10 --use-alleles-trigger false --disable-optimizations false --just-determine-active-regions false --dont-genotype false --dont-trim-active-regions false --max-disc-ar-extension 25 --max-gga-ar-extension 300 --padding-around-indels 150 --padding-around-snps 20 --kmer-size 10 --kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --recover-dangling-heads false --do-not-recover-dangling-branches false --min-dangling-branch-length 4 --consensus false --max-num-haplotypes-in-population 128 --error-correct-kmers false --min-pruning 2 --debug-graph-transformations false --kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --likelihood-calculation-engine PairHMM --base-quality-score-threshold 18 --pair-hmm-gap-continuation-penalty 10 --pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --native-pair-hmm-threads 4 --native-pair-hmm-use-double-precision false --debug false --use-filtered-reads-for-annotations false --bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --capture-assembly-failure-bam false --error-correct-reads false --do-not-run-physical-phasing false --min-base-quality-score 10 --smith-waterman JAVA --use-new-qual-calculator false --annotate-with-num-discovered-alleles false --heterozygosity 0.001 --indel-heterozygosity 1.25E-4 --heterozygosity-stdev 0.01 --standard-min-confidence-threshold-for-calling 10.0 --max-genotype-count 1024 --sample-ploidy 2 --genotyping-mode DISCOVERY --contamination-fraction-to-filter 0.0 --output-mode EMIT_VARIANTS_ONLY --all-site-pls false --min-assembly-region-size 50 --max-assembly-region-size 300 --assembly-region-padding 100 --max-reads-per-alignment-start 50 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --interval-set-rule UNION --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --disable-tool-default-read-filters false --minimum-mapping-quality 20",Version=4.0.2.1,Date="March 26, 2019 3:59:59 PM AWST">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the negative training set of bad variants">
##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the positive training set of good variants">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds of being a true variant versus being false under the trained gaussian mixture model">
##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr4,length=191154276>
##contig=<ID=chr5,length=180915260>
##contig=<ID=chr6,length=171115067>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chr8,length=146364022>
##contig=<ID=chr9,length=141213431>
##contig=<ID=chr10,length=135534747>
##contig=<ID=chr11,length=135006516>
##contig=<ID=chr12,length=133851895>
##contig=<ID=chr13,length=115169878>
##contig=<ID=chr14,length=107349540>
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr16,length=90354753>
##contig=<ID=chr17,length=81195210>
##contig=<ID=chr18,length=78077248>
##contig=<ID=chr20,length=63025520>
##contig=<ID=chrY,length=59373566>
##contig=<ID=chr19,length=59128983>
##contig=<ID=chr22,length=51304566>
##contig=<ID=chr21,length=48129895>
##contig=<ID=chr6_ssto_hap7,length=4928567>
##contig=<ID=chr6_mcf_hap5,length=4833398>
##contig=<ID=chr6_cox_hap2,length=4795371>
##contig=<ID=chr6_mann_hap4,length=4683263>
##contig=<ID=chr6_apd_hap1,length=4622290>
##contig=<ID=chr6_qbl_hap6,length=4611984>
##contig=<ID=chr6_dbb_hap3,length=4610396>
##contig=<ID=chr17_ctg5_hap1,length=1680828>
##contig=<ID=chr4_ctg9_hap1,length=590426>
##contig=<ID=chr1_gl000192_random,length=547496>
##contig=<ID=chrUn_gl000225,length=211173>
##contig=<ID=chr4_gl000194_random,length=191469>
##contig=<ID=chr4_gl000193_random,length=189789>
##contig=<ID=chr9_gl000200_random,length=187035>
##contig=<ID=chrUn_gl000222,length=186861>
##contig=<ID=chrUn_gl000212,length=186858>
##contig=<ID=chr7_gl000195_random,length=182896>
##contig=<ID=chrUn_gl000223,length=180455>
##contig=<ID=chrUn_gl000224,length=179693>
##contig=<ID=chrUn_gl000219,length=179198>
##contig=<ID=chr17_gl000205_random,length=174588>
##contig=<ID=chrUn_gl000215,length=172545>
##contig=<ID=chrUn_gl000216,length=172294>
##contig=<ID=chrUn_gl000217,length=172149>
##contig=<ID=chr9_gl000199_random,length=169874>
##contig=<ID=chrUn_gl000211,length=166566>
##contig=<ID=chrUn_gl000213,length=164239>
##contig=<ID=chrUn_gl000220,length=161802>
##contig=<ID=chrUn_gl000218,length=161147>
##contig=<ID=chr19_gl000209_random,length=159169>
##contig=<ID=chrUn_gl000221,length=155397>
##contig=<ID=chrUn_gl000214,length=137718>
##contig=<ID=chrUn_gl000228,length=129120>
##contig=<ID=chrUn_gl000227,length=128374>
##contig=<ID=chr1_gl000191_random,length=106433>
##contig=<ID=chr19_gl000208_random,length=92689>
##contig=<ID=chr9_gl000198_random,length=90085>
##contig=<ID=chr17_gl000204_random,length=81310>
##contig=<ID=chrUn_gl000233,length=45941>
##contig=<ID=chrUn_gl000237,length=45867>
##contig=<ID=chrUn_gl000230,length=43691>
##contig=<ID=chrUn_gl000242,length=43523>
##contig=<ID=chrUn_gl000243,length=43341>
##contig=<ID=chrUn_gl000241,length=42152>
##contig=<ID=chrUn_gl000236,length=41934>
##contig=<ID=chrUn_gl000240,length=41933>
##contig=<ID=chr17_gl000206_random,length=41001>
##contig=<ID=chrUn_gl000232,length=40652>
##contig=<ID=chrUn_gl000234,length=40531>
##contig=<ID=chr11_gl000202_random,length=40103>
##contig=<ID=chrUn_gl000238,length=39939>
##contig=<ID=chrUn_gl000244,length=39929>
##contig=<ID=chrUn_gl000248,length=39786>
##contig=<ID=chr8_gl000196_random,length=38914>
##contig=<ID=chrUn_gl000249,length=38502>
##contig=<ID=chrUn_gl000246,length=38154>
##contig=<ID=chr17_gl000203_random,length=37498>
##contig=<ID=chr8_gl000197_random,length=37175>
##contig=<ID=chrUn_gl000245,length=36651>
##contig=<ID=chrUn_gl000247,length=36422>
##contig=<ID=chr9_gl000201_random,length=36148>
##contig=<ID=chrUn_gl000235,length=34474>
##contig=<ID=chrUn_gl000239,length=33824>
##contig=<ID=chr21_gl000210_random,length=27682>
##contig=<ID=chrUn_gl000231,length=27386>
##contig=<ID=chrUn_gl000229,length=19913>
##contig=<ID=chrM,length=16571>
##contig=<ID=chrUn_gl000226,length=15008>
##contig=<ID=chr18_gl000207_random,length=4262>
##source=ApplyVQSR
##source=CombineGVCFs
##source=GenotypeGVCFs
##source=HaplotypeCaller
##bcftools_isecVersion=1.7+htslib-1.7
##bcftools_isecCommand=isec -p intersect ./variants/sam_merged.vcf.gz ./variants/genotype_all.snp_recalibrated.indel_recalibrated.bgzip.vcf.gz; Date=Fri Mar 29 07:12:42 2019
##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">
##INFO=<ID=OLD_VARIANT,Number=.,Type=String,Description="Original chr:pos:ref:alt encoding">
##VEP="v97" time="2019-11-15 10:57:53" cache="/data/workspace/richard/software/vep_cache/homo_sapiens/97_GRCh37" ensembl=97.378db18 ensembl-io=97.dc917e1 ensembl-variation=97.26a059c ensembl-funcgen=97.24f4d3c 1000genomes="phase3" COSMIC="86" ClinVar="201810" ESP="20141103" HGMD-PUBLIC="20174" assembly="GRCh37.p13" dbSNP="151" gencode="GENCODE 19" genebuild="2011-04" gnomAD="r2.1" polyphen="2.2.2" regbuild="1.0" sift="sift5.2.2"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|HGVS_OFFSET">
##bcftools_viewVersion=1.9-222-g7171da0+htslib-1.9-293-g5e83884-dirty
##bcftools_viewCommand=view -f PASS combined_HaplotypeCaller.d.n.vep.vcf; Date=Fri Nov 15 11:43:44 2019
##dbNSFPCommand="java -Xmx5g search_dbNSFP40a -i /data/workspace/richard/zika/tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf -o /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out -v hg19 > /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out.log 2>&1"
##INFO=<ID=ACMG,Number=.,Type=String,Description="ACMG classification and evidence. Format: classification|pvs1|ps1|pm1|pm2|pm4|pm5|pp2|pp3">
##bcftools_annotateVersion=1.9-222-g7171da0+htslib-1.9-293-g5e83884-dirty
##bcftools_annotateCommand=annotate -a combined_dbnsfp_geminiPASSonly_acmg_anno.tsv.gz -h combined_dbnsfp_geminiPASSonly_acmg_anno.hdr.tsv -c CHROM,POS,REF,ALT,ACMG -o ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf -O v ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf; Date=Mon Nov 18 10:54:25 2019
##bcftools_viewCommand=view -S ^samples2exclude.txt -O z -o combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg_45samples.vcf.gz combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf.gz; Date=Sat Jan 18 20:51:50 2020
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	JN-102	JN-104	JN-11	JN-110	JN-114	JN-118	JN-12	JN-120	JN-122	JN-126	JN-128	JN-13	JN-136	JN-138DIL	JN-146	JN-149	JN-154DIL	JN-163	JN-165	JN-167	JN-170	JN-172	JN-174	JN-176	JN-187	JN-198DIL2	JN-206	JN-213	JN-217	JN-235	JN-239	JN-252LAB	JN-262	JN-264	JN-267	JN-283	JN-31	JN-329	JN-34	JN-344	JN-52	JN-71	JN-91	JN-94	JN-96DIL
```
