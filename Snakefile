configfile: "config.json"

rule merge_lanes:
	input:
		expand("output/{lane}.sorted.dedup.recal.realigned.bam", lane=config["lanes"].keys())
	output:
		"output/{samplename}_merged.bam"
	shell:
		"echo 'picard MergeBamFiles'"

rule bwa_map_lane:
	input:
		config["reference"],
		lambda wildcards: config["lanes"][wildcards.lane]
	output:
		"output/{lane}.bam"
	shell:
		"bwa mem {input} | samtools view -Sb - > {output}"

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
		"output/{lane}.sorted.dedup.bam"
	output:
		"output/{lane}.sorted.dedup.recal.bam"
	shell:
		"echo 'recalibrating'"

rule gatk_realign_indels:
	input:
		"output/{lane}.sorted.dedup.recal.bam"
	output:
		"output/{lane}.sorted.dedup.recal.realigned.bam"
	shell:
		"echo realigning indels"


