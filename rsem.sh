#input star.transcriptome_bam

for bam in /processed_files/*Aligned.toTranscriptome.out.bam;
do
	echo $bam
	prefix="/processed_files/"
	suffix=".Aligned.toTranscriptome.out.bam"
	string=${bam#$prefix}
	string=${string%$suffix}
	echo $string
	rsem-calculate-expression --num-threads 8 --fragment-length-max 1000 --no-bam-output --paired-end --strandedness reverse --estimate-rspd --bam $bam /resources/hg38_rnaseq_refs/rsem_reference_GRCh38_gencode26_ercc/rsem_reference ${string}.rsem
	echo "done!"
done
