## zika_processing
# important note that vt has been altered in ~/software/vt
# see https://github.com/atks/vt/issues/46 for the change required

# original vcf data obtained from Alexia. This had been run through GATK and samtools and the product is the intersect from those tools. This data is in ./raw
# this vcf file and an ANNOVAR annotated version are in the folder raw_all_vcf
# there are 660,007 variants in these files.
```{bash}
ls raw
```

# variant data had not been decomposed or normalised so this was performed via the seqnextgen pipeline, see sng_script.sh
```{bash}
cat sng_script.sh
```

# this copies the original vcf into ./tmp and normalises and decomposes with vt
# the normalised VCF is then annotated with VEP
# the final file combined_HaplotypeCaller.d.n.vep.vcf.gz is a bgzipped, tabix and grabix indexed VCF containing 661866 variants (1859 more than the original due to the normalisation process)
# note that the data is loaded to GEMINI with the --passonly flag and so is equivalent to the combined_HaplotypeCaller.d.n.vep.PASSonly.vcf file compared to dbNSFP below
```{bash}
ls tmp
```

# filter the VCF to PASS only
# this file conains variants that were deemed false positives by GATK (see VQSR docs https://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr)
# it was agreed that this data should only contain PASSed variants so filter this data accordingly and index
# combined_HaplotypeCaller.d.n.vep.PASSonly.vcf contains 612761 variants
bcftools view -f PASS combined_HaplotypeCaller.d.n.vep.vcf > combined_HaplotypeCaller.d.n.vep.PASSonly.vcf
bgzip combined_HaplotypeCaller.d.n.vep.PASSonly.vcf 
tabix combined_HaplotypeCaller.d.n.vep.PASSonly.vcf.gz
grabix index combined_HaplotypeCaller.d.n.vep.PASSonly.vcf.gz

# compare this file to dbNSFP to find ns and sv variants
# note that it is not clear whether or not dbNSFP is normalised and decomposed however doing this with a non-normalised file matched 110 fewer variants
# the dbNSFP data is installed in /data/workspace/alexia/dbNSFP and the following was executed in this directory to allow access ot the dbNSFP suite 
java -Xmx5g search_dbNSFP40a -i /data/workspace/richard/zika/tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf -o /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out -v hg19 > /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out.log 2>&1

# the output of the above produced
# 58111 SNP(s) are found. Written to /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out
# 554650 SNP(s) are not found. Written to /data/workspace/richard/zika/zika_vep_dbnsfp.PASSonly.out.err
# however the above PASSonly.out file contains 58556 lines of variant data (445 more than stated due to duplicates (e.g. 19:44919189 and 16:1992649) arising from alternative consequences in different transcripts)

# to get the required annotation for input into the acmg tool the VCF was brought into GEMINI
# relevant data fields were manually extracted from within the resulting sqlite database
sqlite3 Illumina.db
.headers on
.mode tabs
.output zika_gemini_export.tsv
select chrom,start,end,ref,alt,clinvar_sig,impact,gene,exon,vep_hgvsc,vep_hgvsp,aa_change,codon_change,pfam_domain,max_aaf_all,polyphen_score,sift_score,cadd_scaled from variants;

# find exported GEMINI data that was found in dbNSFP
# use GEMINI chr start + 1 == sbNSFP hg19_pos(1-based)
# dbnsfp has 470 columns so just narrow it down to the hg19 columns we want. include the ref and alt to make sure we are matching the right variant found in dbNSFP.
# "ref" "alt" "hg19_chr" "hg19_pos(1-based)"
```{bash}
awk 'BEGIN{FS="\t"; OFS="\t";}{print $3,$4,$8,$9; }' zika_vep_dbnsfp.PASSonly.out > zika_vep_dbnsfp.PASSonly.hg19cols.out
```

# use R to do the mergin and output a file for input into acmg
# this file contains 58556 lines of variant data representing the subset of the GEMINI data that also matched to dbNSFP.
# note that this includes the 445 duplicates mentioned above
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

# run ACMG
# there are 58556 lines in the input file which will take a while to run
# as there is interaction with live databases across networks, things can go wrong
# to save re-running all 60K variants the file was split into 30 files containing roughly 2K variants each and a job script made
# once completed combine the files
# these files are in ./splitfiles
# note there were 6 variants that had connection errors during the process. Their data (see linesWithConnectionErrors.txt) were manually replaced in all_acmg.out.tsv folloing re-running with acmg
```{bash}
split -l 2000 combined_dbnsfp_geminiPASSonly.tsv acmg_job_
for i in acmg*; do sed -i '1s/^/chrom\tstart\tnewstart\tend\tref\talt\tclinvar_sig\timpact\tgene\texon\tvep_hgvsc\tvep_hgvsp\taa_change\tcodon_change\tpfam_domain\tmax_aaf_all\tpolyphen_score\tsift_score\tcadd_scaled\n/' $i; done
for i in acmg*; do echo "python3 ../../software/ppTimo/Phenoparser/scripts/acmg.py -s $i -p $i -g /data/scratch/reference/acmg/acmg_id.db > $i"_acmg.log" 2>&1"; done > runacmg.sh
chmod +x runacmg.sh
./runacmg.sh
head -n 1 acmg_job_aaacmg.out.tsv > combined_dbnsfp_geminiPASSonly_acmg.tsv
awk 'BEGIN{FS="\t";} FNR>1 {print;}' acmg_job_*acmg.out.tsv >> combined_dbnsfp_geminiPASSonly_acmg.tsv
```

# make a new VCF file of this subset of variants to annotate
# first use zika_vep_dbnsfp.PASSonly.hg19cols.out to make a file to search for these variants in the original VCF
# then grep out the header lines and dbNSFP variants
# note that this file contains 58111 variants which excludes the 445 duplicate findings in dbNSFP mentioned above
```{bash}
awk 'BEGIN{FS="\t";OFS="\t";} FNR>1 {print "chr"$3,$4,".",$1,$2; }' zika_vep_dbnsfp.PASSonly.hg19cols.out > getdbNSFPsubset.txt
grep -e "^#" tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf > ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf
grep -f getdbNSFPsubset.txt tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf >> ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf
```

# edit ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf to add the dbNSFP search command
```{bash}
```

# make a file containing the acmg annotation for these variants to add to the VCF annotation
```{bash}
awk 'BEGIN{FS="\t";OFS="\t";} FNR >1 {print "chr"$1,$3,$5,$6,$20"|"$22"|"$24"|"$26"|"$28"|"$30"|"$32"|"$34"|"$36; }' combined_dbnsfp_geminiPASSonly_acmg.tsv | sort -k1,1V -k2,2n > combined_dbnsfp_geminiPASSonly_acmg_anno.tsv
bgzip combined_dbnsfp_geminiPASSonly_acmg_anno.tsv
tabix -s1 -b2 -e2 combined_dbnsfp_geminiPASSonly_acmg_anno.tsv.gz
```

# annotate dbNSFP variant VCF with ACMG
```{bash}
echo '##INFO=<ID=ACMG,Number=.,Type=String,Description="ACMG classification and evidence. Format: classification|pvs1|ps1|pm1|pm2|pm4|pm5|pp2|pp3">' > combined_dbnsfp_geminiPASSonly_acmg_anno.hdr.tsv
bcftools annotate -a combined_dbnsfp_geminiPASSonly_acmg_anno.tsv.gz -h combined_dbnsfp_geminiPASSonly_acmg_anno.hdr.tsv -c CHROM,POS,REF,ALT,ACMG -o ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf -O v ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.vcf
bgzip ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf
tabix ./tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.dbNSFPonly.acmg.vcf.gz
```
