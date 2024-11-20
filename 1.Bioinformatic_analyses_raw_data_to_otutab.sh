#!/bin/bash
#BSUB -q Q96C1T_X12
#BSUB -R span[hosts=1]
#BSUB -J cl_job.sh
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 10

echo `hostname`
# date
echo job start time is `date`
cd ~/Eufungi/

# Raw reads data in ~/Eufungi/rawreadsfile/

# step0: quality-checked，检查序列质量
mkdir -p rawreadsfile/fastqc/
# Filenames=($(find rawreadsfile/ -type f -wholename "*.fq.gz" -exec bash -c 'echo $(basename {} .fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
# for i in $(seq 0 $((${#Filenames[@]} - 1))); do
# echo ${Filenames[i]}
# done
# zcat rawreadsfile/${Filenames[i]}_R1.fq.gz | head
# zcat rawreadsfile/${Filenames[i]}_R2.fq.gz | head

Filenames=($(find rawreadsfile/ -type f -wholename "*.fq.gz" -exec bash -c 'echo $(basename {} .fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
for i in $(seq 0 $((${#Filenames[@]} - 1))); do
echo ${Filenames[i]}
fastqc rawreadsfile/${Filenames[i]}_R1.fq.gz -o rawreadsfile/fastqc/ -t 10
fastqc rawreadsfile/${Filenames[i]}_R2.fq.gz -o rawreadsfile/fastqc/ -t 10
done



# step1: demultiplexed，样品拆分
mkdir -p step1_demul/
Filenames=($(find rawreadsfile/ -type f -wholename "*.fq.gz" -exec bash -c 'echo $(basename {} .fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
for i in $(seq 0 $((${#Filenames[@]} - 1))); do
fastq-multx -B barcodes.txt -m 0 -b rawreadsfile/${Filenames[i]}_R1.fq.gz rawreadsfile/${Filenames[i]}_R2.fq.gz -o step1_demul/${Filenames[i]}_%.R1.fq.gz -o step1_demul/${Filenames[i]}_%.R2.fq.gz -q 20 -x
done
rm step1_demul/*_unmatched.R*.fq.gz



# step2: remove low-quality 3'-ends，fastp 3'端 质控
mkdir -p step2_qualitycontrol/
mkdir -p step2_qualitycontrol/fastplog/
Filenames=($(find step1_demul/ -type f -name "*R1.fq.gz" -exec bash -c 'echo $(basename {} .R1.fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
for i in $(seq 0 `expr ${#Filenames[@]} - 1`); do
fastp -i step1_demul/${Filenames[i]}.R1.fq.gz -I step1_demul/${Filenames[i]}.R2.fq.gz \
  -c -A -f 0 -F 0 -q 15 -u 40 -n 5 -3 -W 4 --cut_tail_window_size 4 --cut_tail_mean_quality 20 -l 100 --thread 10 \
  -o step2_qualitycontrol/${Filenames[i]}.R1.fq.gz -O step2_qualitycontrol/${Filenames[i]}.R2.fq.gz \
  -h step2_qualitycontrol/fastplog/${Filenames[i]}_fastppreport.html -j step2_qualitycontrol/fastplog/${Filenames[i]}_fastppreport.json -R step2_qualitycontrol/fastplog/${Filenames[i]}_fastppreport
done
# zcat step1_demul/${Filenames[i]}.R1.fq.gz | head
# zcat step1_demul/${Filenames[i]}.R2.fq.gz | head
# zcat step2_qualitycontrol/${Filenames[i]}.R1.fq.gz | head
# zcat step2_qualitycontrol/${Filenames[i]}.R2.fq.gz | head

# mkdir -p step2_qualitycontrol/fastqc/
# Filenames=($(find step2_qualitycontrol/ -type f -name "*R1.fq.gz" -exec bash -c 'echo $(basename {} .R1.fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
# for i in $(seq 0 `expr ${#Filenames[@]} - 1`); do
# fastqc step2_qualitycontrol/${Filenames[i]}.R1.fq.gz -o step2_qualitycontrol/fastqc/ -t 10
# fastqc step2_qualitycontrol/${Filenames[i]}.R2.fq.gz -o step2_qualitycontrol/fastqc/ -t 10
# done



# step3: 双端序列合并，将文件压缩，大幅度节省内存空间
mkdir -p step3_merged/
Filenames=($(find step2_qualitycontrol/ -type f -name "*R1.fq.gz" -exec bash -c 'echo $(basename {} .R1.fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
for i in $(seq 0 `expr ${#Filenames[@]} - 1`); do
vsearch --fastq_mergepairs step2_qualitycontrol/${Filenames[i]}.R1.fq.gz --reverse step2_qualitycontrol/${Filenames[i]}.R2.fq.gz \
  --threads 10 --relabel ${Filenames[i]}. --fastqout - \
  | gzip -6 > step3_merged/${Filenames[i]}.fq.gz
done
zcat step3_merged/${Filenames[i]}.fq.gz | head
seqkit stat step3_merged/*.fq.gz > step3_merged/seqkit_all.txt

# mkdir -p step3_merged/fastqc/
# Filenames=($(find step3_merged/ -type f -name "*.fq.gz" -exec bash -c 'echo $(basename {} .fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
# for i in $(seq 0 `expr ${#Filenames[@]} - 1`); do
# fastqc step3_merged/${Filenames[i]}.fq.gz -o step3_merged/fastqc/ -t 10
# done

Filenames=($(find step3_merged/ -type f -name "*.fq.gz" -exec bash -c 'echo $(basename {} .fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
for i in $(seq 0 `expr ${#Filenames[@]} - 1`); do
zcat step3_merged/${Filenames[i]}.fq.gz >> step3_merged/FIPMe_all.fq
done
gzip -6 step3_merged/FIPMe_all.fq
# mkdir -p step3_merged/fastqc/
# fastqc step3_merged/FIPMe_all.fq.gz -o step3_merged/fastqc/ -t 10
seqkit stat step3_merged/FIPMe_all.fq.gz >> step3_merged/seqkit_all.txt

# 删除临时文件
find step3_merged/ -type f -name "*.fq.gz" ! -name "FIPMe_all.fq.gz" -exec rm -f {} +
rm -r step1_demul/
rm -r step2_qualitycontrol/



# step4: trimming Barcodes and primers，样品序列切除引物与质控，barcodes和引物切除，--fastq_maxee_rate 0.01 质控
# R: A, G; Y: C, T
# GTGARTCATCGAATCTTTG...TCCTCCGCTTATTGATATGC
# GTGAATCATCRAATYTTTG...TCCTCCGCTTATTGATATGC
# GTGARTCATCRAATYTTTG...GCATATCAATAAGCGGAGGA
mkdir -p step4_cutadapt/
mkdir -p step4_cutadapt/untrim/
# Filenames=($(find step3_merged/ -type f -name "*.fq.gz" -exec bash -c 'echo $(basename {} .fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
# for i in $(seq 0 `expr ${#Filenames[@]} - 1`); do
# cutadapt step3_merged/${Filenames[i]}.fq.gz \
  # -g "GTGARTCATCRAATYTTTG;max_error_rate=0.25;required"..."GCATATCAATAAGCGGAGGA;max_error_rate=0.25;required" \
  # --output step4_cutadapt/${Filenames[i]}_trim.fq.gz \
  # --untrimmed-output step4_cutadapt/untrim/${Filenames[i]}_untrim.fq.gz
# done
# zcat step4_cutadapt/${Filenames[i]}_trim.fq.gz | head
# seqkit stat step4_cutadapt/*_trim.fq.gz > step4_cutadapt/seqkit_all_trim.txt

cutadapt step3_merged/FIPMe_all.fq.gz \
  -g "GTGARTCATCRAATYTTTG;max_error_rate=0.25;required"..."GCATATCAATAAGCGGAGGA;max_error_rate=0.25;required" \
  --output step4_cutadapt/FIPMe_all_trim.fq.gz \
  --untrimmed-output step4_cutadapt/untrim/FIPMe_all_untrim.fq.gz
seqkit stat step4_cutadapt/*_trim.fq.gz >> step4_cutadapt/seqkit_all_trim.txt



# step5: reads with ≥1 ambiguous nucleotide or >1% expected error rate，质控
mkdir -p step5_filtered/
# Filenames=($(find step4_cutadapt/ -type f -name "*_trim.fq.gz" -exec bash -c 'echo $(basename {} _trim.fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
# for i in $(seq 0 `expr ${#Filenames[@]} - 1`); do
# vsearch --fastx_filter step4_cutadapt/${Filenames[i]}_trim.fq.gz \
  # --fastq_maxns 0 --fastq_maxee_rate 0.01 \
  # --fastq_minlen 100 --fastqout - \
  # | gzip -6 > step5_filtered/${Filenames[i]}_filtered.fq.gz
# done
# zcat step5_filtered/${Filenames[i]}_filtered.fq.gz | head
# seqkit stat step5_filtered/*_filtered.fq.gz > step5_filtered/seqkit_all_filtered.txt

# Filenames=($(find step4_cutadapt/ -type f -name "*_trim.fq.gz" -exec bash -c 'echo $(basename {} _trim.fq.gz)' \; | sed 's/_R[12]$//' | awk '!seen[$0]++'))
# for i in $(seq 0 `expr ${#Filenames[@]} - 1`); do
# zcat step5_filtered/${Filenames[i]}_filtered.fq.gz >> step5_filtered/FIPMe_all_filtered.fq
# done
# gzip -6 step5_filtered/FIPMe_all_filtered.fq
# seqkit stat step5_filtered/FIPMe_all_filtered.fq.gz > step5_filtered/seqkit_FIPMe_all_filtered.txt

vsearch --fastx_filter step4_cutadapt/FIPMe_all_trim.fq.gz \
  --fastq_maxns 0 --fastq_maxee_rate 0.01 \
  --fastqout - \
  | gzip -6 > step5_filtered/FIPMe_all_filtered.fq.gz
seqkit stat step5_filtered/FIPMe_all_filtered.fq.gz > step5_filtered/seqkit_FIPMe_all_filtered.txt



# step6: extract TS2 region using ITSxpress，提取ITS2序列，去除ITS2序列长度小于100bp的序列，将文件压缩，大幅度节省内存空间
mkdir -p step6_itsxpress/
itsxpress --fastq step5_filtered/FIPMe_all_filtered.fq.gz \
  --cluster_id 1.0 --single_end --region "ITS2" --taxa "All" \
  --outfile step6_itsxpress/FIPMe_all_filtered_ITS2.fq.gz \
  --log step6_itsxpress/FIPMe_all_filtered_ITS2.log \
  --threads 10
seqkit stat step6_itsxpress/FIPMe_all_filtered_ITS2.fq.gz > step6_itsxpress/seqkit_FIPMe_all_filtered_ITS2.txt

# vsearch --fastx_filter step6_itsxpress/FIPMe_all_filtered_ITS2.fq.gz \
  # --fastq_minlen 100 --fasta_width 0 --threads 10 \
  # --fastaout step6_itsxpress/FIPMe_all_filtered_ITS2.fa

vsearch --fastx_filter step6_itsxpress/FIPMe_all_filtered_ITS2.fq.gz \
  --fastq_minlen 100 --fasta_width 0 --threads 10 \
  --fastaout - \
  | gzip -6 > step6_itsxpress/FIPMe_all_filtered_ITS2.fa.gz
seqkit stat step6_itsxpress/FIPMe_all_filtered_ITS2.fa.gz > step6_itsxpress/seqkit_FIPMe_all_filtered_ITS2_fa.txt



# # 提取Uchime ITS2数据库，准备用于基于参考数据库的嵌合体去除（可选，Uchime有已经提取好的ITS2数据库）
# ITSx -i unite_dataset/uchime_reference_dataset_16_10_2022.fasta \
  # --complement T --save_regions all --graphical F --positions T \
  # -E 1e-1 -t all --cpu 10 --preserve T \
  # -o unite_dataset/uchime_reference_dataset_16_10_2022

# # 提取Unite（和RDP ITS) ITS2数据库，准备用于分类分配的数据库
# ITSx -i unite_dataset/unite_dataset_all_dataset_25.07.2023.fa \
  # --complement T --save_regions all --graphical F --positions T \
  # -E 1e-1 -t all --cpu 10 --preserve T \
  # -o unite_dataset/unite_dataset_all_dataset_25.07.2023



# step7: 序列去冗余
# 参数--miniuniqusize设置使用序列的最小出现频次，默认为8，此处设置为10。
mkdir -p step7_uniques/
# vsearch --derep_fulllength step6_itsxpress/FIPMe_all_filtered_ITS2.fa \
  # --minuniquesize 10 --sizeout --relabel Uni_ \
  # --fasta_width 0 --threads 10 \
  # --output step7_uniques/EufungiMe_all_uniques.fa

vsearch --derep_fulllength step6_itsxpress/FIPMe_all_filtered_ITS2.fa.gz \
  --minuniquesize 10 --sizeout --relabel Uni_ \
  --fasta_width 0 --threads 10 \
  --output step7_uniques/EufungiMe_all_uniques.fa



# step8: uniques去嵌合体（先基于数据库去嵌合、再自身比对去嵌合更好），以备聚类OTU
mkdir -p step8_dchime/
vsearch --uchime_ref step7_uniques/EufungiMe_all_uniques.fa \
  --db unite_dataset/uchime_reference_dataset_16_10_2022/uchime_reference_dataset_16_20_2022_ITS2.fasta \
  --fasta_width 0 --sizein --sizeout --threads 10 \
  --uchimeout step8_dchime/EufungiMe_all_uniques_ref.txt --uchimeout5 \
  --chimeras step8_dchime/EufungiMe_all_uniques_chimeras_ref.fa \
  --nonchimeras step8_dchime/EufungiMe_all_uniques_nonchimeras_ref.fa

vsearch --uchime3_denovo step8_dchime/EufungiMe_all_uniques_nonchimeras_ref.fa \
  --fasta_width 0 --sizein --sizeout --threads 10 \
  --uchimeout step8_dchime/EufungiMe_all_uniques_ref_denovo.txt --uchimeout5 \
  --chimeras step8_dchime/EufungiMe_all_uniques_chimeras_ref_denovo.fa \
  --nonchimeras step8_dchime/EufungiMe_all_uniques_nonchimeras_ref_denovo.fa



# step9: 聚类OTU
mkdir -p step9_otus/
vsearch --cluster_size step8_dchime/EufungiMe_all_uniques_nonchimeras_ref_denovo.fa  \
  --id 0.97 --relabel OTU_ --fasta_width 0 \
  --qmask none --sizein --sizeout \
  --centroids step9_otus/EufungiMe_all_otus_size.fa

# 删除序列频率
# sed -i 's/;.*//' 
sed 's/;.*//' step9_otus/EufungiMe_all_otus_size.fa > step9_otus/EufungiMe_all_otus.fa



# step10: 生成OTU特征表
mkdir -p step10_otustabs/
# vsearch --usearch_global step6_itsxpress/FIPMe_all_filtered_ITS2.fa \
  # --db step9_otus/EufungiMe_all_otus.fa --id 0.97 --threads 10 \
  # --otutabout step10_otustabs/EufungiMe_all_otutab.txt

vsearch --usearch_global step6_itsxpress/FIPMe_all_filtered_ITS2.fa.gz \
  --db step9_otus/EufungiMe_all_otus.fa --id 0.97 --threads 10 \
  --otutabout step10_otustabs/EufungiMe_all_otutab.txt



# step11: OTU sintax
mkdir -p step11_otusintax/
vsearch --sintax step9_otus/EufungiMe_all_otus.fa \
  --db unite_dataset/unite_dataset_all_dataset_25.07.2023/ITSX.ITS2.fasta \
  --sintax_cutoff 0.8 \
  --tabbedout step11_otusintax/EufungiMe_all_otus_unite.sintax



# step12: OTU unite taxonomy
mkdir -p step12_taxonomy/
cut -f 1,4 step11_otusintax/EufungiMe_all_otus_unite.sintax \
  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
  > step12_taxonomy/EufungiMe_all_otus_unite_taxonomy2.txt

awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
  step12_taxonomy/EufungiMe_all_otus_unite_taxonomy2.txt > step12_taxonomy/EufungiMe_all_otus_unite_otus.tax

sed 's/;/\t/g;s/.__//g;' step12_taxonomy/EufungiMe_all_otus_unite_otus.tax|cut -f 1-8 | \
  sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
  > step12_taxonomy/EufungiMe_all_otus_unite_taxonomy.txt



# step11: OTU blastn
blastn -query step9_otus/EufungiMe_all_otus.fa -db ~/blast_db/NT/nt/nt -evalue 1e-30 -task blastn \
  -best_hit_score_edge 0.05 -best_hit_overhang 0.25 \
  -max_target_seqs 10 -max_hsps 1 -num_threads 10 \
  -outfmt "6 qaccver qlen sgi saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms nident qcovs qcovhsp qcovus" \
  -out step11_otusintax/EufungiMe_all_otus_nt.txt



# step12: OTU blast taxonomy
# 去重，仅保留第一个blastn结果
awk '!seen[$1]++' step11_otusintax/EufungiMe_all_otus_nt.txt > step12_taxonomy/EufungiMe_all_otus_nt_1.txt

cut -f 1,15-24 step12_taxonomy/EufungiMe_all_otus_nt_1.txt > step12_taxonomy/EufungiMe_all_otus_nt_1.txt_taxid.txt

for i in `cat step12_taxonomy/EufungiMe_all_otus_nt_1.txt_taxid.txt | csvtk -t cut -f 2`; do
echo ${i} | taxonkit lineage -j 4 | taxonkit reformat -P -F -t -j 4 | cut -f 1,3 >> step12_taxonomy/EufungiMe_all_otus_nt_1.txt_taxid_taxon.txt
done

paste -d '\t' step12_taxonomy/EufungiMe_all_otus_nt_1.txt_taxid.txt  step12_taxonomy/EufungiMe_all_otus_nt_1.txt_taxid_taxon.txt > step12_taxonomy/EufungiMe_all_otus_nt_1.txt_taxid_taxon_lineage.txt


# 去重，仅保留第二个blastn结果
awk 'seen[$1]++ == 1' step11_otusintax/EufungiMe_all_otus_nt.txt > step12_taxonomy/EufungiMe_all_otus_nt_2.txt

cut -f 1,15-24 step12_taxonomy/EufungiMe_all_otus_nt_2.txt > step12_taxonomy/EufungiMe_all_otus_nt_2.txt_taxid.txt

for i in `cat step12_taxonomy/EufungiMe_all_otus_nt_2.txt_taxid.txt | csvtk -t cut -f 2`; do
echo ${i} | taxonkit lineage -j 4 | taxonkit reformat -P -F -t -j 4 | cut -f 1,3 >> step12_taxonomy/EufungiMe_all_otus_nt_2.txt_taxid_taxon.txt
done

paste -d '\t' step12_taxonomy/EufungiMe_all_otus_nt_2.txt_taxid.txt  step12_taxonomy/EufungiMe_all_otus_nt_2.txt_taxid_taxon.txt > step12_taxonomy/EufungiMe_all_otus_nt_2.txt_taxid_taxon_lineage.txt


# 去重，仅保留第三个blastn结果
awk 'seen[$1]++ == 2' step11_otusintax/EufungiMe_all_otus_nt.txt > step12_taxonomy/EufungiMe_all_otus_nt_3.txt

cut -f 1,15-24 step12_taxonomy/EufungiMe_all_otus_nt_3.txt > step12_taxonomy/EufungiMe_all_otus_nt_3.txt_taxid.txt

for i in `cat step12_taxonomy/EufungiMe_all_otus_nt_3.txt_taxid.txt | csvtk -t cut -f 2`; do
echo ${i} | taxonkit lineage -j 4 | taxonkit reformat -P -F -t -j 4 | cut -f 1,3 >> step12_taxonomy/EufungiMe_all_otus_nt_3.txt_taxid_taxon.txt
done

paste -d '\t' step12_taxonomy/EufungiMe_all_otus_nt_3.txt_taxid.txt  step12_taxonomy/EufungiMe_all_otus_nt_3.txt_taxid_taxon.txt > step12_taxonomy/EufungiMe_all_otus_nt_3.txt_taxid_taxon_lineage.txt



# date
echo job end time is `date`
echo this job is cl_job.sh



# Reference:
# Andrews S. FastQC: a quality control tool for high throughput sequence data. https://github.com/s-andrews/FastQC (2010).
# Aronesty E. Comparison of Sequencing Utility Programs. The Open Bioinformatics Journal 7, 1-8 (2013).
# Rognes T, Flouri T, Nichols B, Quince C, Mahe F. VSEARCH: a versatile open source tool for metagenomics. Peerj 4, e2584 (2016).
# Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884-i890 (2018).
# Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet journal 17, 10-12 (2011).
# Edgar RC, Flyvbjerg H. Error filtering, pair assembly and error correction for next-generation sequencing reads. Bioinformatics 31, 3476-3482 (2015).
# Rivers A, Weber K, Gardner T, Liu S, Armstrong S. ITSxpress: Software to rapidly trim internally transcribed spacer sequences with quality scores for marker gene analysis. https://f1000research.com/articles/7-1418/v1 (2018). https://github.com/USDA-ARS-GBRU/itsxpress
# Nilsson RH, et al. An open source software package for automated extraction of ITS1 and ITS2 from fungal ITS sequences for use in high-throughput community assays and molecular ecology. Fungal Ecol 4,284–287 (2010).
# Bengtsson-Palme J, et al. Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods Ecol Evol 10,914–919 (2013).
# Abarenkov K, et al. UNITE UCHIME reference dataset. Version 9.0.). UNITE Community (2022).
# Tedersoo L, et al. Best practices in metabarcoding of fungi: From experimental design to results. Mol Ecol 31, 2769-2795 (2022).
# Abarenkov K, et al. UNITE USEARCH/UTAX release for eukaryotes. Version 9.0.). UNITE Community (2023).
# Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome 6, 226 (2018).
# Nguyen NH, Smith D, Peay K, Kennedy P. Parsing ecological signal from noise in next generation amplicon sequencing. New Phytol 205, 1389-1393 (2015).
# Cobian GM, Egan CP, Amend AS. Plant–microbe specificity varies as a function of elevation. The ISME Journal 13, 2778-2788 (2019).
# The Global Soil Mycobiome consortium dataset: https://github.com/Mycology-Microbiology-Center/GSMc
# Easy Amplicon data analysis pipeline: https://github.com/YongxinLiu/EasyAmplicon
# Pipelines for analyzing metabarcoding data: https://github.com/BPerezLamarque/Scripts
# R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing. https://www.R-project.org/ (2013).



################################################################################
################################################################################
################################################################################