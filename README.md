# MixA-seq
1. Create an empty folder
2. Prepare fast_5_barcode
 >fast_barcode
 GACTGCGTACG
3. Prepare fast_3_barcode
 >p3-K3-33
 TGGTGTCAAACC
 >p3-K17-39
 TGGATACAAACC
 >p3-sample name
 TATTATCAAACC
#-Sequence name must start with `p3-`
#-Followed by a combination of characters (sample name), the characters only include `-`, `_`, numbers, capital letters
4. Prepare the reference genome fasta, fai, dict
#-bwa index ref.fasta
#-picard CreateSequenceDictionary R=ref.fasta
5. "ln -s /local_path/bin .; cp /local_path/run.sh ." to the current project folder
6. Modify the parameters in "run.sh" and execute "sh run.sh" when everything is ready
Note: A "block" is separated by a blank line. If you want to skip a certain block, you can use "#" to comment out the entire block, but the first "block" is the basic configuration and must be retained.
Result: Statistical file "plot.tsv" for the number of SNPs of all samples
Result: Single sample SNP chromosome mapping "sample name.diff.png" and "sample name.diff.pdf", pdf is a vector diagram, which can be easily modified in Adobe Illustrator
