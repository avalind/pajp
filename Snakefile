import os

configfile: "config.json"

def extract_picard_readgroup(identifier):
	parts = identifier.split("_")
	dataset = {"lane": parts[0], "date": parts[1], "flowcell": parts[2], "scilife_id": parts[3]+"_"+parts[4]}
	id_ = dataset["flowcell"] + "_" + dataset["lane"]
	return "RGID=%s RGSM=%s RGPL=illumina RGLB=%s RGPU=%s" % \
		(id_, dataset["scilife_id"], dataset["scilife_id"], id_)

rule merge_lanes:
	input:
		expand("output/{file_id}.sorted.dedup.realigned.recal.bam", file_id=config["lanes"].keys())
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

rule picard_add_or_replace_read_groups:
	input:
		"output/{file}.bam"
	output:
		"output/{file}.sorted.bam"
	params:
		java_args_picard="",
		readgroup_picard=lambda wildcards: extract_picard_readgroup(wildcards.file)
	shell:
		"java {params.java_args_picard} -jar $PICARD_HOME/AddOrReplaceReadGroups.jar "
		"I={input} "
		"O={output} "
		"SO=coordinate {params.readgroup_picard} "
		
rule picard_dedup:
	input:
		"output/{file}.sorted.bam"
	output:
		"output/{file}.sorted.dedup.bam",
		"metadata/{file}.dedup_metrics.txt"
	params:
		java_args_picard="",
	shell:
		"java {params.java_args_picard} -jar $PICARD_HOME/MarkDuplicates.jar "
		"I={input} O={output[0]} METRICS_FILE={output[1]} "

rule gatk_realign_indels_make_targets:
	input:
		"output/{file}.sorted.dedup.bam"
	output:
		"metadata/{file}.target_intervals.list"
	params:
		gold_indel_mills="",
		gold_indel_1000g="",
		gatk_refpath="",
		java_args_gatk="",
	shell:
		"java {params.java_args_gatk} -jar $GATK_HOME/GenomeAnalysisTK.jar "
		"-T RealignerTargetCreator "
		"-R {params.gatk_refpath} "
		"-I {input} "
		"-known {params.gold_indel_mills} -known {params.gold_indel_1000g} "
		"-o {output}"

rule gatk_realign_indels_apply_targets:
	input:
		"output/{file}.sorted.dedup.bam",
		"metadata/{file}.target_intervals.list"
	output:
		"output/{file}.sorted.dedup.realigned.bam"
	params:
		gold_indel_mills="",
		gold_indel_1000g="",
		gatk_refpath="",
		java_args_gatk="",
	shell:
		"java {params.java_args_gatk} -jar $GATK_HOME/GenomeAnalysisTK.jar"
		"-T IndelRealigner "
		"-R {params.gatk_refpath} "
		"-I {input[0]} "
		"-targetIntervals {input[1]} "
		"-known {params.gold_indel_mills} -known {params.gold_indel_1000g} "
		"-o {output} "		

rule gatk_recalibrate_bsqr:
	input:
		"output/{file}.sorted.dedup.realigned.bam"
	output:
		"output/{file}.sorted.dedup.realigned.recal.bam"
	shell:
		"echo recalculate"


