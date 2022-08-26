
## Refine the reference genome for 5' CAR-T alignment, including Kyrmiah and Yescarta

# Inputs: 
ref=hg38_gencode34_CAR_v0 # Directory with cellranger reference genome containing CAR guesses
rundir="CAR_assembly/" # Directory with fastqs containing CAR sequence

cd $rundir

# STEP 1 - Align to approximate CAR reference
STAR --genomeDir $ref/star \
	--runThreadN 16  --readFilesIn $(ls -1 *.fastq.gz | tr '\n' ',') \
	--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat
samtools index Aligned.sortedByCoord.out.bam

## STEP 2 - Call differences in actual sequencing from guess
samtools view -h Aligned.sortedByCoord.out.bam Yescarta Kymriah | \
	awk -F '\t' 'BEGIN{OFS=FS}{gsub("N","D",$6);print}' | \
	bcftools mpileup -f $ref/fasta/genome.fa /dev/stdin  | \
	bcftools call -mv | \
	bcftools filter -e 'IMF < .5' | \
	bgzip > plasmid_variants.vcf.gz
tabix plasmid_variants.vcf.gz

## STEP 3 - Create Consensus fasta and updated gtf
cat $ref/fasta/genome.fa | \
	bcftools consensus --chain liftover.chain plasmid_variants.vcf.gz > reference_updated.fa
CrossMap.py gff liftover.chain hg38_gencode34_CAR_v0/genes/genes.gtf genes_updated.gtf

## STEP 4 - Re-build STAR index
cellranger mkref --genome=hg38_gencode34_CAR5p_v1 \
	        --fasta=reference_updated.fa \
		        --genes=genes_updated.gtf \
			        --nthreads=16 \
				        --memgb=100
## Zip up final reference
tar czvf -I pigz hg38_gencode34_CAR5p_v1.tar.gz hg38_gencode34_CAR5p_v1


