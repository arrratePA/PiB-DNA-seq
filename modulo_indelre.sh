#!/bin/bash
############################################################################################################################
## DNAseq_Fastq_to_annotated variants
###########################################################################################################################
############################################################################################################################
#script configuration
    #colors
BOLDGREEN="\e[1;32m"
ITALICSRED="\e[3;31m"
ENDCOLOR="\e[0m"
    # exit when any command fails
set -e
    # keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG # no funciona
    # echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT #no funciona

#0. Define aliases
alias picard="java -jar /home/bioaraba/bioinfo/bioinfo_tools/picard.jar"
alias gatk3="java -jar /home/bioaraba/bioinfo/bioinfo_tools/gatk3/GenomeAnalysisTK.jar"
#0. Define default variables
humandb="/home/bioaraba/bioinfo/bioinfo_tools/annovar/humandb/"
dict="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"
dbsnp138="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
hg38_indels="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
gs_indels="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

#0. Get variables defined by input
while getopts w:g:i:s: flag
do
    case "${flag}" in
        w) workingpath=$(realpath ${OPTARG});;
        g) genome_path=$(realpath ${OPTARG});;
        i) intervals=$(realpath ${OPTARG});;
        s) samples=$(realpath ${OPTARG});;
    esac
done
echo "Working path: $workingpath";
echo "genome_path: $genome_path";
echo "intervals: $intervals";
echo "samples: $samples";

#0. Create list of samples and analysis folders
cd $workingpath
mkdir bam_files_indelre
mkdir vcf_files_indelre
mkdir annotation_files_indelre
bam_files_indelre=$(realpath bam_files_indelre)
vcf_files_indelre=$(realpath vcf_files_indelre)
annotation_files_indelre=$(realpath annotation_files_indelre)


for i in `cat $samples`; do

# Indelrealigner
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}IndelRealigner in ${i} sample${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
    gatk3 -T RealignerTargetCreator \
        -ip 150 \
        -R $genome_path \
        -L $intervals \
        -I $workingpath/bam_files/${i}.RG.sorted.mkdup.bam \
        -o ${i}.RG.sorted.mkdup.bam.realignertargetcreator.intervals

    gatk3 \
        -T IndelRealigner \
        -ip 150 \
        -R $genome_path \
        -targetIntervals ${i}.RG.sorted.mkdup.bam.realignertargetcreator.intervals \
        -I $workingpath/bam_files/${i}.RG.sorted.mkdup.bam \
        -o ${i}.RG.sorted.mkdup.indelrealigned.bam

#5. Base Quality Recalibration (BQSR) indelrealigned bams
    gatk BaseRecalibrator \
        -I ${i}.RG.sorted.mkdup.indelrealigned.bam \
        -R $genome_path \
        -L $intervals \
        -ip 150 \
        -O ${i}.RG.sorted.mkdup.indelrealigned.bqsr.table \
        --known-sites $dbsnp138 \
        --known-sites $hg38_indels \
        --known-sites $gs_indels

    gatk ApplyBQSR \
        -I ${i}.RG.sorted.mkdup.indelrealigned.bam \
        -R $genome_path \
        -L $intervals \
        -ip 150 \
        --bqsr-recal-file ${i}.RG.sorted.mkdup.indelrealigned.bqsr.table \
        -O ${i}.RG.sorted.mkdup.indelrealigned.bqsr.bam

#6. HaplotypeCaller
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Doing HaplotypeCaller in realigned ${i} sample${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log    
    gatk HaplotypeCaller \
        -R $genome_path \
        -I ${i}.RG.sorted.mkdup.indelrealigned.bqsr.bam \
        -L $intervals \
        -ip 150 \
        -O ${i}.indelrealigned.g.vcf.gz \
        -ERC GVCF \
        -bamout ${i}.RG.sorted.mkdup.indelrealigned.bqsr.HC.bamout.bam


    mv *.bam $bam_files_indelre
    mv *.bai $bam_files_indelre
    mv ${i}.indelrealigned.g.vcf* $vcf_files_indelre
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variant Calling of idelrealigment ${i} sample finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

done

#7. Consolidate GVCFs: GenomicsDBImport
#This step consists of consolidating the contents of GVCF files across multiple samples in order to improve scalability and speed the next step, joint genotyping. Note that this is NOT equivalent to the joint genotyping step; variants in the resulting merged GVCF cannot be considered to have been called jointly.
#7.1. Create cohort.sample_map
cd $vcf_files_indelre

for i in *.vcf.gz; do echo `bcftools query -l $i`;echo $i;done | paste - - > cohort.sample_map
#7.1. GenomicsDBImport
gatk GenomicsDBImport \
    --genomicsdb-workspace-path database  \
    --batch-size 8 \
    -L $intervals \
    -ip 150 \
    --sample-name-map cohort.sample_map \
    --reader-threads 4

#7.2. Joint-Call Cohort: GenotypeGVCFs
#At this step, we gather all the per-sample GVCFs (or combined GVCFs if we are working with large numbers of samples) and pass them all together to the joint genotyping tool, GenotypeGVCFs. 
#This produces a set of joint-called SNP and indel calls ready for filtering. This cohort-wide analysis empowers sensitive detection of variants even at difficult sites, and produces a squared-off matrix of genotypes that provides information about all sites of interest in all samples considered, which is important for many downstream analyses.
gatk GenotypeGVCFs \
    -R $genome_path \
    -V gendb://database \
    -O cohort.vcf.gz
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variant Calling of cohort indelrealigment finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
gunzip cohort.vcf.gz

#8. Spliting by sample and filtering by: quality, coverage and 0/0 genotypes
#8.1 Split vcf by sample
for i in `cat $samples`; do
    cd $vcf_files_indelre
    gatk SelectVariants \
        -V cohort.vcf \
        -R $genome_path \
        -sn ${i} \
        --keep-original-ac \
        -O ${i}.vcf

#8.2 Filtering and Quality control
    bcftools +fill-tags ${i}.vcf -o ${i}.tag.vcf  #add new tags, as VAF
    vcftools --vcf ${i}.tag.vcf --minQ 30 --min-meanDP 10 --minGQ 20 --recode --recode-INFO-all --out ${i}_F1 #filter for variants will at least Q30, DP20 and GQ 20
    bcftools view -e 'GT="RR"' ${i}_F1.recode.vcf -o ${i}_filtered.vcf #exclude homozygous for reference
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variants filtering of ${i} sample finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

#8.3 Annotation of SNVs and indels    
    table_annovar.pl ${i}_filtered.vcf $humandb -buildver hg38 -out ${i}_SNVs.myanno -remove --intronhgvs 40 -protocol refGene,cytoBand,clinvar_20150629,revel,dbnsfp42c,dbnsfp30a,dbnsfp31a_interpro,gnomad211_exome,gnomad211_genome,exac03,dbscsnv11,avsnp142,regsnpintron -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
    mv ${i}_SNVs.myanno.hg38_multianno.txt ${i}_SNVs.myanno.hg38_multianno.vcf $annotation_files_indelre
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variants annotation of ${i} sample finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

#8.4 Final table construction
    cd $annotation_files_indelre
    awk '{print $1 ":" $2 ":" $3 ":" $4 ":"  $5 "\t" $0}' ${i}_SNVs.myanno.hg38_multianno.txt > ${i}_SNVs.myanno.hg38_multianno_0.txt
    sed -i "s/^/${i}\t/" ${i}_SNVs.myanno.hg38_multianno_0.txt
    awk '{print $(NF)}' ${i}_SNVs.myanno.hg38_multianno_0.txt | awk -F '[:]' '{print $1, $2, $3, $(NF) }' OFS='\t' | sed '1d' | sed '1i GT\tAD\tDP\tVAF1' > ${i}_SNVs.myanno.hg38_multianno_format.txt
    awk '{print $(NF-2)}' ${i}_SNVs.myanno.hg38_multianno_0.txt | awk -F '[;]' '{print $2, $4, $6}' OFS="\t" | sed 's+=+    +g' | sed '1d' | awk '{print $2, $4, $6}' OFS="\t" | sed '1i AC_Cohort\tAF_Cohort\tAN_Cohort' > ${i}_SNVs.myanno.hg38_multianno_info.txt
    cut -f 1-7 ${i}_SNVs.myanno.hg38_multianno_0.txt > ${i}_final_table1.txt
    cut -f 8- ${i}_SNVs.myanno.hg38_multianno_0.txt > ${i}_final_table4.txt
    paste ${i}_final_table1.txt ${i}_SNVs.myanno.hg38_multianno_format.txt ${i}_SNVs.myanno.hg38_multianno_info.txt ${i}_final_table4.txt > ${i}_indelre_variants_table.txt
    cat ${i}_indelre_variants_table.txt > ${i}_indelre_variants_table.xls
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variants table of ${i} sample created${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
done

#Create file with results for all samples
ls *_indelre_variants_table.txt > samples_SNVs_Library.txt
vim -c "%s/_variants_table.txt//g|wq" samples_SNVs_Library.txt
for i in $(cat samples_SNVs_Library.txt); do
    cat ${i}_variants_table.txt | sed '1d' >> Library_SNVs_multianno.txt
done
sed -n 1p ${i}_variants_table.txt > header
cut -f 2- header | awk '{print $1="Sample\t" $0}' > header_ok 
cat header Library_SNVs_multianno.txt > Library_SNVs_multianno_header.txt

for i in $(cat samples_SNVs_Library.txt); do
    N=$(sed '1d' ${i}_variants_table.txt | wc -l)
    printf "${i}\t${N}\n" | tee -a SNV_count.txt
done

#Remove not needed files
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Removing not needed files${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
cd $vcf_files_indelre
rm -r *_F1.log
rm -r *_SNVs.myanno.avinput
rm -r *.vcf.idx
rm -r *.vcf.gz.tbi
rm -r database
rm -r cohort.*
cd $bam_files_indelre
rm -r *.RG.sorted.mkdup.indelrealigned.bam
rm -r *.RG.sorted.mkdup.indelrealigned.bai
cd $annotation_files_indelre
rm -r *_SNVs.myanno.hg38_multianno_format.txt
rm -r *_SNVs.myanno.hg38_multianno_info.txt
rm -r *_final_table4.txt
rm -r *_final_table1.txt
rm -r header_ok
cd $workingpath
rm -r *.realignertargetcreator.intervals
rm -r *.indelrealigned.bqsr.table

echo -e `date +[%D-%R]` "\t${BOLDGREEN}DNAseq pipeline with indel realigment finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

