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
adapters="/home/bioaraba/bioinfo/bioinfo_tools/cutadapt/adapters/TruSeq3-PE-2.fa"
bed="/home/bioaraba/bioinfo/bioinfo_tools/panelbeds/iPPSD-v1_TE-95220955_hg38_sorted.bed"
intervals="NULL"
gene_list="/home/bioaraba/bioinfo/bioinfo_tools/panelbeds/iPPSD_v1.genelist.refseq"
humandb="/home/bioaraba/bioinfo/bioinfo_tools/annovar/humandb/"
dict="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"
genome_path="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
dbsnp138="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
hg38_indels="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
gs_indels="/home/bioaraba/bioinfo/bioinfo_tools/genome_references/hg38/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
check_report=false #no funciona, si pongo true no lo escribe en la pregunta como true pero luego si me pregunta 
indelre_do=true
pwd_sh=$(realpath $(dirname $0))

#0. Get variables defined by input
while getopts w:a:b:i:d:g:h:r: flag
do
    case "${flag}" in
        w) workingpath=$(realpath ${OPTARG});;
        a) adapters=$(realpath ${OPTARG});;
        b) bed=$(realpath ${OPTARG});;
        i) intervals=$(realpath ${OPTARG});;
        d) dict=$(realpath ${OPTARG});;
        g) genome_path=$(realpath ${OPTARG});;
        h) humandb=$(realpath ${OPTARG});;
        r) check_report=$(realpath ${OPTARG});;
    esac
done

#0. Advise defined variables and paths
echo "Working path: $workingpath";
echo "adapters: $adapters";
echo "bed: $bed";
echo "intervals: $intervals";
echo "dict: $dict";
echo "genome_path: $genome_path";
echo "humandb: $humandb";
echo "pwd_sh: $pwd_sh";
echo "Do you want to ask for checking fastqc report before continuing with aligment?:" $check_report;
echo "Do you want to try the variant calling with indel realigment step?:" $indelre_do;

#0. Create list of samples, analysis folders and real paths
cd $workingpath
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Starting DNAseq pipeline from fastq to annotated variants${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
ls *_L001_R1_001.fastq.gz > samples.txt
vim -c "%s/_L001_R1_001.fastq.gz//g|wq" samples.txt
samples=$(realpath samples.txt)
mkdir quality_files
mkdir bam_files
mkdir vcf_files
mkdir annotation_files
mkdir coverage_statistics
quality_files=$(realpath quality_files)
bam_files=$(realpath bam_files)
vcf_files=$(realpath vcf_files)
annotation_files=$(realpath annotation_files)
coverage_statistics=$(realpath coverage_statistics)

echo -e `date +[%D-%R]` "\t${BOLDGREEN}Directories created${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

#0. Copy bed and create intervals file if not defined as input
if [ "$intervals" = "NULL" ] ; then
    cp $bed .
    used_bed=$(basename $bed)
    picard BedToIntervalList -I $used_bed -O $used_bed.interval_list -SD $dict
    intervals=$(realpath $used_bed.interval_list)
    echo  "interval_list: $intervals"
fi

#1.A. Check quality
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Starting sample quality check and preprocessing${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
for i in `cat $samples`; do
	echo -e `date +[%D-%R]` "\t${BOLDGREEN}Starting analysis of ${i} sample${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

    fastqc ${i}_L001_R1_001.fastq.gz -o $quality_files
    fastqc ${i}_L001_R2_001.fastq.gz -o $quality_files    

#2. Trimming and dropping reads for adapters, quality and length: -j cores, -m min-length and -q quality 
    cutadapt -j 4 -m 100 -q 30 -b file:$adapters -B file:$adapters  ${i}_L001_R1_001.fastq.gz ${i}_L001_R2_001.fastq.gz -o ${i}.1_trim.fq.gz -p ${i}.2_trim.fq.gz
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Trimming step of ${i} sample finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

#1.B. Check quality again after trimming
    fastqc ${i}.1_trim.fq.gz -o $quality_files
    fastqc ${i}.2_trim.fq.gz -o $quality_files

done

#1.C. Create a merge quality report with fastqc analysis of all samples 
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Creating quality check merged report${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
multiqc $quality_files -o $quality_files -i $(basename $workingpath)
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Quality merged report created${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

if $check_report; then
    #open $quality_files/*_multiqc_report.html
    read -p "\t${ITALICSRED}Press [Enter] if fastqc is OK and you want to continue, otherwise CTRL+C to exit${ENDCOLOR}"
fi

#3. Alignment with reference genome.
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Starting aligment to reference $(basename $genome_path)${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

for i in `cat $samples`; do
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Starting aligment of ${i} sample${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

#3.1 Decompress
    gunzip ${i}.1_trim.fq.gz
    gunzip ${i}.2_trim.fq.gz

#3.2 Bwa mem = the algorithm works by seeding alignments with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW).
    bwa mem -t 4 -MP -Y -R '@RG\tID:${i}\tSM:${i}\tLB:Twist\tPL:ILLUMINA' $genome_path ${i}.1_trim.fq ${i}.2_trim.fq > ${i}.sam
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Alignment of ${i} sample finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
    #Create bam file and sort + index
    samtools view -bS ${i}.sam > ${i}.bam
    samtools sort ${i}.bam -o ${i}.sorted.bam
    #AddOrReplaceReadGroups (Picard) --> Assigns all the reads in a file to a single new read-group
    #--RGLB / -LB	Read-Group library
    #--RGPL / -PL	Read-Group platform (e.g. ILLUMINA, SOLID)
    #--RGPU / -PU	Read-Group platform unit (eg. run barcode)
    #--RGSM / -SM	Read-Group sample name
    #To see RG in bam file --> samtools view -H input.bam
    picard AddOrReplaceReadGroups -I ${i}.sorted.bam -O ${i}.RG.sorted.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM ${i}
    #samtools index ${i}.RG.sorted.bam
    rm ${i}.sam
    rm ${i}.bam   
#4. MarkDuplicates
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Doing MarkDuplicates ${i} sample${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
    picard MarkDuplicates I=${i}.RG.sorted.bam O=${i}.RG.sorted.mkdup.bam M=${i}.mkdup_metrics.txt
#index bam
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Doing indexing ${i} sample${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
    samtools index ${i}.RG.sorted.mkdup.bam

#5. Base Quality Recalibration (BQSR) no realigned bams
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Doing BQSR ${i} sample${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
    gatk BaseRecalibrator \
        -I ${i}.RG.sorted.mkdup.bam \
        -R $genome_path \
        -L $intervals \
        -ip 150 \
        -O ${i}.RG.sorted.mkdup.bqsr.table \
        --known-sites $dbsnp138 \
        --known-sites $hg38_indels \
        --known-sites $gs_indels

    gatk ApplyBQSR \
        -I ${i}.RG.sorted.mkdup.bam \
        -R $genome_path \
        -L $intervals \
        -ip 150 \
        --bqsr-recal-file ${i}.RG.sorted.mkdup.bqsr.table \
        -O ${i}.RG.sorted.mkdup.bqsr.bam        

#6. HaplotypeCaller
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Starting HaplotypeCaller without realigned in ${i} sample${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log    
    gatk HaplotypeCaller \
        -R $genome_path \
        --dbsnp $dbsnp138 \
        -I ${i}.RG.sorted.mkdup.bqsr.bam \
        -L $intervals \
        -ip 150 \
        -O ${i}.g.vcf.gz \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        -bamout ${i}.RG.sorted.mkdup.bqsr.HC.bamout.bam
    
    mv *.bam $bam_files
    mv *.bai $bam_files
    mv ${i}.g.vcf* $vcf_files
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variant Calling of ${i} sample finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
done

#7. Consolidate GVCFs: GenomicsDBImport
#This step consists of consolidating the contents of GVCF files across multiple samples in order to improve scalability and speed the next step, joint genotyping. Note that this is NOT equivalent to the joint genotyping step; variants in the resulting merged GVCF cannot be considered to have been called jointly.
#7.1. Create cohort.sample_map
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Consolidate GVCFs and create cohort.sample_map${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
cd $vcf_files

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
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    -O cohort.vcf.gz
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variant Calling of cohort finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
gunzip cohort.vcf.gz

#8. Spliting by sample and filtering by: quality, coverage and 0/0 genotypes
#8.1 Split vcf by sample
for i in `cat $samples`; do
    cd $vcf_files
    gatk SelectVariants \
        -V cohort.vcf \
        -R $genome_path \
        -sn ${i} \
        --keep-original-ac \
        -select "QUAL > 20.0 && DP >= 10 && GQ > 20 (isHomVar == 1) (isHet == 1)"
        -O ${i}.vcf

#8.2 Filtering VCFs by Quality, depth, and genotype 
    bcftools +fill-tags ${i}.vcf -o ${i}.tag.vcf  #add new tags, as VAF

#8.3 Annotation of SNVs and indels    
    table_annovar.pl ${i}.tag.vcf $humandb -buildver hg38 -out ${i}_SNVs.myanno -remove --intronhgvs 40 -protocol refGene,cytoBand,clinvar_20150629,revel,dbnsfp42c,dbnsfp30a,dbnsfp31a_interpro,gnomad211_exome,gnomad211_genome,exac03,dbscsnv11,avsnp142,regsnpintron -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
    mv ${i}_SNVs.myanno.hg38_multianno.txt ${i}_SNVs.myanno.hg38_multianno.vcf $annotation_files
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variants annotation of ${i} sample finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

#8.4 Final table construction
    cd $annotation_files
    awk '{print $1 ":" $2 ":" $3 ":" $4 ":"  $5 "\t" $0}' ${i}_SNVs.myanno.hg38_multianno.txt > ${i}_SNVs.myanno.hg38_multianno_0.txt
    sed -i "s/^/${i}\t/" ${i}_SNVs.myanno.hg38_multianno_0.txt
    awk '{print $(NF)}' ${i}_SNVs.myanno.hg38_multianno_0.txt | awk -F '[:]' '{print $1, $2, $3, $(NF) }' OFS='\t' | sed '1d' | sed '1i GT\tAD\tDP\tVAF1' > ${i}_SNVs.myanno.hg38_multianno_format.txt
    awk '{print $(NF-2)}' ${i}_SNVs.myanno.hg38_multianno_0.txt | awk -F '[;]' '{print $2, $4, $6}' OFS="\t" | sed 's+=+    +g' | sed '1d' | awk '{print $2, $4, $6}' OFS="\t" | sed '1i AC_Cohort\tAF_Cohort\tAN_Cohort' > ${i}_SNVs.myanno.hg38_multianno_info.txt
    cut -f 1-7 ${i}_SNVs.myanno.hg38_multianno_0.txt > ${i}_final_table1.txt
    cut -f 8- ${i}_SNVs.myanno.hg38_multianno_0.txt > ${i}_final_table4.txt
    paste ${i}_final_table1.txt ${i}_SNVs.myanno.hg38_multianno_format.txt ${i}_SNVs.myanno.hg38_multianno_info.txt ${i}_final_table4.txt > ${i}_variants_table.txt
    cat ${i}_variants_table.txt > ${i}_variants_table.xls
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Variants table of ${i} sample created${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
done

#Create file with results for all samples
ls *_variants_table.txt > samples_SNVs_Library.txt
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

#9. Calculate coverage statistics 
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Calculating coverage statistics${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
cd $bam_files
for i in `cat $samples`; do
    gatk DepthOfCoverage \
        -R $genome_path \
        -O $coverage_statistics/${i} \
        -I ${i}.RG.sorted.mkdup.bam \
        -L $intervals \
        -ip 50 \
        -gene-list $gene_list \
        --summary-coverage-threshold 20 \
        --summary-coverage-threshold 50 \
        --summary-coverage-threshold 100
done

echo -e `date +[%D-%R]` "\t${BOLDGREEN}DNAseq pipeline finished${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log

#Remove not needed files
echo -e `date +[%D-%R]` "\t${BOLDGREEN}Removing not needed files${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
cd $vcf_files
rm -r *_F1.log
rm -r *_SNVs.myanno.avinput
rm -r *.vcf.idx
rm -r *.vcf.gz.tbi
rm -r database
cd $bam_files
rm -r *.sorted.bam 
cd $quality_files
rm -r *_fastqc.zip
cd $annotation_files
rm -r *_SNVs.myanno.hg38_multianno_format.txt
rm -r *_SNVs.myanno.hg38_multianno_info.txt
rm -r *_final_table4.txt
rm -r *_final_table1.txt
rm -r header
cd $workingpath
rm -r *_trim.fq
rm -r *_metrics.txt
rm -r *.bqsr.table


if $indelre_do; then
    echo -e `date +[%D-%R]` "\t${BOLDGREEN}Doing variant calling with indel realigment step${ENDCOLOR}" | tee -a $workingpath/timeElapseReport.log
    cd $workingpath
    sh $pwd_sh/modulo_indelre.sh -w $workingpath -i $intervals -g $genome_path -s $samples
fi
