##1. Create an empty folder
##2. Prepare fast_5_barcode
## >fast_barcode
## GACTGCGTACG
##3. Prepare fast_3_barcode
## >p3-K3-33
## TGGTGTCAAACC
## >p3-K17-39
## TGGATACAAACC
## >p3-sample name
## TATTATCAAACC
##-Sequence name must start with `p3-`
##-Followed by a combination of characters (sample name), the characters only include `-`, `_`, numbers, capital letters
##4. Prepare the reference genome fasta, fai, dict
##-bwa index ref.fasta
##-picard CreateSequenceDictionary R=ref.fasta
##5. "ln -s /local_path/bin .; cp /local_path/run.sh ." to the current project folder
##6. Modify the parameters in "run.sh" and execute "sh run.sh" when everything is ready
##Note: A "block" is separated by a blank line. If you want to skip a certain block, you can use "#" to comment out the entire block, but the first "block" is the basic configuration and must be retained.
##Result: Statistical file "plot.tsv" for the number of SNPs of all samples
##Result: Single sample SNP chromosome mapping "sample name.diff.png" and "sample name.diff.pdf", pdf is a vector diagram, which can be easily modified in Adobe Illustrator
RAW1="clean/customer-0lszfITV/mushuhunhe19-5_FDDP190650516-1a/mushuhunhe19-5_FDDP190650516-1a_1.clean.fq.gz"
RAW2="clean/customer-0lszfITV/mushuhunhe19-5_FDDP190650516-1a/mushuhunhe19-5_FDDP190650516-1a_2.clean.fq.gz"
REFG="03.reference/Mesculenta_305_v6.fasta"
PARENT="K3-1 K7-21" # Separated by spaces，Parents ID；
CHR="Chromosome01,Chromosome02,Chromosome03,Chromosome04,Chromosome05,Chromosome06,Chromosome07,Chromosome08,Chromosome09,Chromosome10,Chromosome11,Chromosome12,Chromosome13,Chromosome14,Chromosome15,Chromosome16,Chromosome17,Chromosome18" # 以逗号分隔，与REFG中的序列名对应
WINDOW=1000 # KB
DRAW=100 # KB
source `which env_parallel.sh`

fastqsplitter -i $RAW1 -o split.01.R1.fq.gz -o split.02.R1.fq.gz -o split.03.R1.fq.gz -o split.04.R1.fq.gz -o split.05.R1.fq.gz -o split.06.R1.fq.gz -o split.07.R1.fq.gz -o split.08.R1.fq.gz -o split.09.R1.fq.gz -o split.10.R1.fq.gz -o split.11.R1.fq.gz -o split.12.R1.fq.gz
fastqsplitter -i $RAW2 -o split.01.R2.fq.gz -o split.02.R2.fq.gz -o split.03.R2.fq.gz -o split.04.R2.fq.gz -o split.05.R2.fq.gz -o split.06.R2.fq.gz -o split.07.R2.fq.gz -o split.08.R2.fq.gz -o split.09.R2.fq.gz -o split.10.R2.fq.gz -o split.11.R2.fq.gz -o split.12.R2.fq.gz

rm -rf tmp tmp3; mkdir -p tmp tmp3
ls split.*.R1.fq.gz | env_parallel --env PATH -j 12 './bin/scanAP -i {} -a fast_5_barcode -s ./tmp/{= s#.fq.gz## =}.stat -d ./tmp/{= s#.fq.gz## =}.detail -k 10 -e 2 && ./bin/trim_seq.pl {} -prefix ./tmp/{= s#.fq.gz## =} -trim_detail ./tmp/{= s#.fq.gz## =}.detail --edge 20 -trim_mode both --len_p 0.3 --len_t 25'
ls split.*.R2.fq.gz | env_parallel --env PATH -j 12 './bin/scanAP -i {} -a fast_5_barcode -s ./tmp/{= s#.fq.gz## =}.stat -d ./tmp/{= s#.fq.gz## =}.detail -k 10 -e 2 && ./bin/trim_seq.pl {} -prefix ./tmp/{= s#.fq.gz## =} -trim_detail ./tmp/{= s#.fq.gz## =}.detail --edge 20 -trim_mode both --len_p 0.3 --len_t 25'
ls split.*.R1.fq.gz | env_parallel --env PATH -j 12 './bin/scanAP -i ./tmp/{= s#.fq.gz## =}.trim.fastq -a fast_3_barcode -s ./tmp3/{= s#.fq.gz## =}.stat -d ./tmp3/{= s#.fq.gz## =}.detail -k 10 -e 2 && ./bin/trim_seq.pl ./tmp/{= s#.fq.gz## =}.trim.fastq -prefix ./tmp3/{= s#.fq.gz## =} -trim_detail ./tmp3/{= s#.fq.gz## =}.detail --edge 20 -trim_mode both --len_p 0.3 --len_t 25'
ls split.*.R2.fq.gz | env_parallel --env PATH -j 12 './bin/scanAP -i ./tmp/{= s#.fq.gz## =}.trim.fastq -a fast_3_barcode -s ./tmp3/{= s#.fq.gz## =}.stat -d ./tmp3/{= s#.fq.gz## =}.detail -k 10 -e 2 && ./bin/trim_seq.pl ./tmp/{= s#.fq.gz## =}.trim.fastq -prefix ./tmp3/{= s#.fq.gz## =} -trim_detail ./tmp3/{= s#.fq.gz## =}.detail --edge 20 -trim_mode both --len_p 0.3 --len_t 25'
ls split.*.R1.fq.gz | env_parallel --env PATH -j 12 './bin/fltfastq2pe.pl -fastq1 ./tmp3/{= s#.R1.fq.gz## =}.R1.trim.fastq -fastq2 ./tmp3/{= s#.R1.fq.gz## =}.R2.trim.fastq'
ls split.*.R1.fq.gz | env_parallel --env PATH -j 12 './bin/proc_head.sh ./tmp3/{= s#.R1.fq.gz## =}.R1.trim.pair.fastq > ./tmp3/{= s#.R1.fq.gz## =}.R1.trim.pair.h.fastq'
ls split.*.R1.fq.gz | env_parallel --env PATH -j 12 './bin/proc_head.sh ./tmp3/{= s#.R1.fq.gz## =}.R2.trim.pair.fastq > ./tmp3/{= s#.R1.fq.gz## =}.R2.trim.pair.h.fastq'
cat tmp3/split.*.R1.trim.pair.h.fastq > all.R1.trim.pair.h.fastq
cat tmp3/split.*.R2.trim.pair.h.fastq > all.R2.trim.pair.h.fastq
awk 'NR%4{printf "%s\t",$0;next;}1' all.R1.trim.pair.h.fastq > all.R1.trim.pair.h.tab
awk 'NR%4{printf "%s\t",$0;next;}1' all.R2.trim.pair.h.fastq > all.R2.trim.pair.h.tab
paste all.R1.trim.pair.h.tab all.R2.trim.pair.h.tab | sed 's# @#\t@#g' | awk -v OFS="\t" -v FS="\t" '{if ($1=="@NONE") print $6,$2,$3,$4,$5,$7,$8,$9,$10; else print $1,$2,$3,$4,$5,$7,$8,$9,$10}' > all.trim.pair.h.tab
rm -rf tmp4; mkdir -p tmp4
awk -v OFS="\t" -v FS="\t" '{print>"./tmp4/"$1}' all.trim.pair.h.tab
ls tmp4/@p3-* | sed 's#tmp4/@p3-##' > sample.list

mkdir -p 01_splitdata
cat sample.list | env_parallel --env PATH -j 12 'cut -f 2-5 ./tmp4/@p3-{} | sed "s#\t#\n#g" | gzip -n > 01_splitdata/{}.R1.fq.gz'
cat sample.list | env_parallel --env PATH -j 12 'cut -f 6-9 ./tmp4/@p3-{} | sed "s#\t#\n#g" | gzip -n > 01_splitdata/{}.R2.fq.gz'
mkdir -p tmp 02_alignment 
conda deactivate
source activate bwa
cat sample.list | env_parallel --env PATH --env REFG -j 12 'bwa mem -t 4 -R "@RG\tID:{}\tSM:{}" $REFG 01_splitdata/{}.R1.fq.gz 01_splitdata/{}.R2.fq.gz | picard SortSam TMP_DIR=./tmp CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate INPUT=/dev/stdin OUTPUT=02_alignment/{}.sorted.bam' &> log.bwa
cat sample.list | env_parallel --env PATH -j 12 'picard MarkDuplicates TMP_DIR=./tmp ASSUME_SORT_ORDER=coordinate INPUT=02_alignment/{}.sorted.bam OUTPUT=02_alignment/{}.sorted.md.bam METRICS_FILE=02_alignment/{}.sorted.md.metrics' &> log.picard
mkdir -p 03_vcf
cat sample.list | env_parallel --env PATH --env REFG -j 12 'gatk HaplotypeCaller -ERC GVCF -R $REFG -I 02_alignment/{}.sorted.bam -O 03_vcf/{}.sorted.vcf.gz' &> log.gatk

conda deactivate
source activate vcftools
(zcat $(ls 03_vcf/*.sorted.vcf.gz | head -1) | head -n 3000 | grep "^#" | sed '$ d' && echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") > header
cat sample.list | env_parallel --env PATH -j 12 'vcftools --gzvcf 03_vcf/{}.sorted.vcf.gz --minDP 3 --minQ 30 --recode --recode-INFO-all --out {}' &> log.filter
ls *recode.vcf | env_parallel --env PATH -j 12 'bgzip -f {}; tabix -f -p vcf {}.gz' && rm *.vcf

##
conda deactivate
source activate vcftools
rename p_recode recode *.p_recode.vcf.gz*
for F in $PARENT; do mv $F.recode.vcf.gz $F.p_recode.vcf.gz; mv $F.recode.vcf.gz.tbi $F.p_recode.vcf.gz.tbi; done
vcftools --gzvcf `echo $PARENT | cut -d " " -f 1`.p_recode.vcf.gz --gzdiff `echo $PARENT | cut -d " " -f 2`.p_recode.vcf.gz --diff-site --out tmp &> log.parent
awk '$4 == "B"' tmp.diff.sites_in_files > tmp.both
vcftools --gzvcf `echo $PARENT | cut -d " " -f 1`.p_recode.vcf.gz --positions-overlap tmp.both --recode --out tmp.both &>> log.parent
mv tmp.both.recode.vcf parent.common.vcf; rm tmp.diff.sites_in_files
for F in $PARENT; do
    ls *.recode.vcf.gz | env_parallel --env PATH --env F -j 12 'vcftools --gzvcf {} --gzdiff $F.p_recode.vcf.gz   --diff-site --out `echo {} | cut -d "." -f 1`:$F' &> log.diff.$F
done
ls *.diff.sites_in_files | env_parallel --env PATH -j 12 'awk '"'"'$4 == "B"'"'"' `echo {} | cut -d "." -f 1`.diff.sites_in_files > `echo {} | cut -d "." -f 1`.both; vcftools --gzvcf `echo {} | cut -d ":" -f 1`.recode.vcf.gz --positions-overlap `echo {} | cut -d "." -f 1`.both --recode --out `echo {} | cut -d "." -f 1`.both; vcftools --gzvcf `echo {} | cut -d "." -f 1`.both.recode.vcf --gzdiff parent.common.vcf  --diff-site --out `echo {} | cut -d "." -f 1`.p; awk '"'"'$4 == "1"'"'"' `echo {} | cut -d "." -f 1`.p.diff.sites_in_files; awk '"'"'$4 == "1"'"'"' `echo {} | cut -d "." -f 1`.p.diff.sites_in_files > `echo {} | cut -d "." -f 1`.p.1; vcftools --gzvcf `echo {} | cut -d ":" -f 1`.recode.vcf.gz --positions-overlap `echo {} | cut -d "." -f 1`.p.1 --recode --out `echo {} | cut -d "." -f 1`.p.1' &> log.vcf
ls *.recode.vcf.gz | env_parallel --env PATH --env PARENT -j 12 'vcftools --gzvcf `echo {} | cut -d "." -f 1`:`echo $PARENT | cut -d " " -f 1`.p.1.recode.vcf --gzdiff `echo {} | cut -d "." -f 1`:`echo $PARENT | cut -d " " -f 2`.p.1.recode.vcf   --diff-site --out `echo {} | cut -d "." -f 1`.pcmp; awk '"'"'$4 == "B"'"'"' `echo {} | cut -d "." -f 1`.pcmp.diff.sites_in_files > `echo {} | cut -d "." -f 1`.pcmp.h; vcftools --gzvcf `echo {} | cut -d "." -f 1`:`echo $PARENT | cut -d " " -f 1`.p.1.recode.vcf --positions-overlap `echo {} | cut -d "." -f 1`.pcmp.h --recode --out `echo {} | cut -d "." -f 1`.pcmp.h; awk '"'"'$4 == "1"'"'"' `echo {} | cut -d "." -f 1`.pcmp.diff.sites_in_files > `echo {} | cut -d "." -f 1`.pcmp.1; vcftools --gzvcf `echo {} | cut -d "." -f 1`:`echo $PARENT | cut -d " " -f 1`.p.1.recode.vcf --positions-overlap `echo {} | cut -d "." -f 1`.pcmp.1 --recode --out `echo {} | cut -d "." -f 1`.pcmp.1; awk -v OFS="\t" '"'"'$4 == "2" {print $1,$3,$2,$4,$6,$5,$8,$7}'"'"' `echo {} | cut -d "." -f 1`.pcmp.diff.sites_in_files > `echo {} | cut -d "." -f 1`.pcmp.2; vcftools --gzvcf `echo {} | cut -d "." -f 1`:`echo $PARENT | cut -d " " -f 2`.p.1.recode.vcf --positions-overlap `echo {} | cut -d "." -f 1`.pcmp.2 --recode --out `echo {} | cut -d "." -f 1`.pcmp.2' &> log.pcmp
rename .pcmp.1.recode :`echo $PARENT | cut -d " " -f 1`.diff *.pcmp.1.recode.vcf; rename .pcmp.2.recode :`echo $PARENT | cut -d " " -f 2`.diff *.pcmp.2.recode.vcf; rename .pcmp.h.recode :hybrid.diff *.pcmp.h.recode.vcf
cut -f 1,2 $REFG.fai > tmp.window && bedtools makewindows -g tmp.window -w $(($WINDOW * 1000)) > windows.bed && rm tmp.window
conda deactivate
source activate r-test
PARENT=$PARENT" hybrid"
echo "SN $PARENT" | sed 's/ /\t/g' > plot.tsv
ls *.recode.vcf.gz | env_parallel --env PATH --env REFG --env PARENT --env CHR --env WINDOW --env DRAW -j 12 'Rscript ./bin/vcf2plot.R -s `echo {} | cut -d "." -f 1` -p `echo $PARENT | sed "s/ /,/g"` -i $REFG.fai -c $CHR -g zq001 -w $WINDOW -d $DRAW' &> log.plot

