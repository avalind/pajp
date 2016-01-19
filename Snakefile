configfile: "config.json"

rule merge_lanes:
	input:
		expand("output/{file_id}.sorted.dedup.recal.realigned.bam", file_id=config["lanes"].keys())
	output:
		"output/{samplename}_merged.bam"
	shell:
		"echo 'picard MergeBamFiles'"

rule bwa_map_lane:
	input:
		reference=config["reference"],
		readpair=lambda wildcards: config["lanes"][wildcards.lane]
	output:
		"output/{lane}.bam"
	shell:
		"""
		trim_galore --fastqc --output_dir trimmed/ --paired {input.readpair} 
		bwa mem {input.reference} trimmed/{wildcards.lane}_1_val_1.fq.gz trimmed/{wildcards.lane}_2_val_2_.fq.gz |
		samtools view -Sb - > {output} 
		"""

rule samtools_sort:
	input:
		"output/{file}.bam"
	output:
		"output/{file}.sorted.bam"
	shell:
		"echo 'samtools sort'"
		
rule picard_dedup:
	input:
		"output/{file}.sorted.bam"
	output:
		"output/{file}.sorted.dedup.bam"
	shell:
		"echo 'deduping'"

rule gatk_recalibrate_bsqr:
	input:
		"output/{file}.sorted.dedup.bam"
	output:
		"output/{file}.sorted.dedup.recal.bam"
	shell:
		"echo 'recalibrating'"

rule gatk_realign_indels:
	input:
		"output/{file}.sorted.dedup.recal.bam"
	output:
		"output/{file}.sorted.dedup.recal.realigned.bam"
	shell:
		"echo realigning indels"


