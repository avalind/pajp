import os

configfile: "test.json"

def all_lanes_for_sample(wildcards):
	skel = "{0}/output/dedup/".format(wildcards.samplename)
	input_files = []
	lane_names = config["samples"][wildcards.samplename].keys()
	for lane_name in lane_names:
		input_files.append(skel+"{0}.sorted.dedup.realigned.recal.bam".format(lane_name))
	return input_files

def get_readpair(wildcards):
	current_sample_name = wildcards.samplename
	current_lane = wildcards.lane
	print(current_lane)
	return config["samples"][current_sample_name][current_lane]	

def extract_picard_readgroup(identifier):
	parts = identifier.split("_")
	dataset = {"lane": parts[0], "date": parts[1], "flowcell": parts[2], "scilife_id": parts[3]+"_"+parts[4]}
	id_ = dataset["flowcell"] + "_" + dataset["lane"]
	return "RGID=%s RGSM=%s RGPL=illumina RGLB=%s RGPU=%s" % \
		(id_, dataset["scilife_id"], dataset["scilife_id"], id_)

rule all:
	input: expand("{samplename}/{samplename}_merged.bam", samplename=config["samples"].keys())
	params:
		time="1:00:00"	

rule merge_lanes:
	input:
		all_lanes_for_sample
	output:
		"{samplename}/{samplename}_merged.bam"
	params:
		java_args_picard=config["java_args_picard"],
		time="8:00:00"
	run:
		input_file_string = ""
		for inp in input:
			input_file_string += " I="+inp
		shell("module load bioinfo-tools java picard && java {params.java_args_picard} -jar $PICARD_HOME/MergeSamFiles.jar O={output} SO=coordinate "+input_file_string)

rule bwa_map_lane:
	input:
		reference=config["bwa_ref"],
		readpair=get_readpair
	output:
		temp("{samplename}/raw_output/{lane}.bam")
	params:
		time="8:00:00"
	log:
		"{samplename}/logs/{lane}.bwa.log"
	shell:
		"""
		module load bioinfo-tools bwa samtools cutadapt FastQC TrimGalore 
		trim_galore --fastqc --output_dir {wildcards.samplename}/trimmed/ --paired {input.readpair} 
		bwa mem -M {input.reference} {wildcards.samplename}/trimmed/{wildcards.lane}_1_val_1.fq.gz {wildcards.samplename}trimmed/{wildcards.lane}_2_val_2.fq.gz 2> {log} |
		samtools view -Sb - > {output} 
		"""

rule picard_add_or_replace_read_groups:
	input:
		"{sample}/raw_output/{file}.bam"
	output:
		temp("{sample}/output/sorted/{file}.sorted.bam")
	params:
		java_args_picard=config["java_args_picard"],
		readgroup_picard=lambda wildcards: extract_picard_readgroup(wildcards.file),
		time="8:00:00"
	shell:
		"""
		module load bioinfo-tools java picard
		java {params.java_args_picard} -jar $PICARD_HOME/AddOrReplaceReadGroups.jar I={input} O={output} SO=coordinate {params.readgroup_picard}
		"""
		
rule picard_dedup:
	input:
		"{sample}/output/sorted/{file}.sorted.bam"
	output:
		temp("{sample}/output/dedup/{file}.sorted.dedup.bam"),
		"{sample}/metadata/{file}.dedup_metrics.txt"
	params:
		java_args_picard=config["java_args_picard"],
		time="8:00:00"
	shell:
		"""
		module load bioinfo-tools java picard
		java {params.java_args_picard} -jar $PICARD_HOME/MarkDuplicates.jar I={input} O={output[0]} METRICS_FILE={output[1]}
		"""

rule samtools_index:
	input:
		"{sample}/output/dedup/{file}.sorted.dedup.bam"
	output:
		"{sample}/output/dedup/{file}.sorted.dedup.bam.bai"
	params:
		time="4:00:00"
	shell:
		"""
		module load bioinfo-tools samtools
		samtools index {input}
		"""

rule gatk_realign_indels_make_targets:
	input:
		"{sample}/output/dedup/{file}.sorted.dedup.bam",
		"{sample}/output/dedup/{file}.sorted.dedup.bam.bai"
	output:
		"{sample}/metadata/{file}.target_intervals.list"
	params:
		gold_indel_mills=config["gold_indel_mills"],
		gold_indel_1000g=config["gold_indel_1000g"],
		gatk_refpath=config["gatk_ref"],
		java_args_gatk=config["java_args_gatk"],
		time="10:00:00"
	shell:
		"""
		module load bioinfo-tools java GATK
		"""
		"java {params.java_args_gatk} -jar $GATK_HOME/GenomeAnalysisTK.jar "
		"-T RealignerTargetCreator "
		"-R {params.gatk_refpath} "
		"-I {input[0]} "
		"-known {params.gold_indel_mills} -known {params.gold_indel_1000g} "
		"-o {output}"

rule gatk_realign_indels_apply_targets:
	input:
		"{sample}/output/dedup/{file}.sorted.dedup.bam",
		"{sample}/metadata/{file}.target_intervals.list"
	output:
		temp("{sample}/output/dedup/{file}.sorted.dedup.realigned.bam")
	params:
		gold_indel_mills=config["gold_indel_mills"],
		gold_indel_1000g=config["gold_indel_1000g"],
		gatk_refpath=config["gatk_ref"],
		java_args_gatk=config["java_args_gatk"],
		time="10:00:00"
	shell:
		"""
		module load bioinfo-tools java GATK 
		"""
		"java {params.java_args_gatk} -jar $GATK_HOME/GenomeAnalysisTK.jar "
		"-T IndelRealigner "
		"-R {params.gatk_refpath} "
		"-I {input[0]} "
		"-targetIntervals {input[1]} "
		"-known {params.gold_indel_mills} -known {params.gold_indel_1000g} "
		"-o {output} "		

rule gatk_recalibrate_calc_bsqr:
	input:
		"{sample}/output/dedup/{file}.sorted.dedup.realigned.bam"
	output:
		"{sample}/metadata/{file}.recal_data.table"
	params:
		gold_indel_mills=config["gold_indel_mills"],
		gold_indel_1000g=config["gold_indel_1000g"],
		dbsnp_138=config["dbsnp_138"],
		gatk_refpath=config["gatk_ref"],
		java_args_gatk=config["java_args_gatk"],
		time="10:00:00"
	shell:
		"""
		module load bioinfo-tools java GATK
		"""
		"java {params.java_args_gatk} -jar $GATK_HOME/GenomeAnalysisTK.jar "
		"-T BaseRecalibrator "
		"-R {params.gatk_refpath} "
		"-I {input} "
		"-knownSites {params.gold_indel_mills} -knownSites {params.gold_indel_1000g} "
		"-knownSites {params.dbsnp_138} "
		"-o {output} "

rule gatk_recalibrate_compare_bsqr:
	input:
		"{sample}/output/dedup/{file}.sorted.dedup.realigned.bam",
		"{sample}/metadata/{file}.recal_data.table",
	output:
		"{sample}/metadata/{file}.post_recal_data.table"
	params:
		gold_indel_mills=config["gold_indel_mills"],
		gold_indel_1000g=config["gold_indel_1000g"],
		dbsnp_138=config["dbsnp_138"],
		gatk_refpath=config["gatk_ref"],
		java_args_gatk=config["java_args_gatk"],
		time="10:00:00"
	shell:
		"""
		module load bioinfo-tools java GATK
		"""
		"java {params.java_args_gatk} -jar $GATK_HOME/GenomeAnalysisTK.jar "
		"-T BaseRecalibrator "
		"-R {params.gatk_refpath} "
		"-I {input[0]} "
		"-knownSites {params.gold_indel_mills} -knownSites {params.gold_indel_1000g} "
		"-knownSites {params.dbsnp_138} "
		"-BQSR {input[1]} "
		"-o {output} "

rule gatk_recalibrate_apply_bsqr:
	input:
		"{sample}/output/dedup/{file}.sorted.dedup.realigned.bam",
		"{sample}/metadata/{file}.recal_data.table",
		"{sample}/metadata/{file}.post_recal_data.table",
	output:
		"{sample}/metadata/{file}.recalibration_plots.pdf",
		temp("{sample}/output/dedup/{file}.sorted.dedup.realigned.recal.bam"),
	params:
		gatk_refpath=config["gatk_ref"],
		java_args_gatk=config["java_args_gatk"],
		time="10:00:00"
	shell:
		"""
		module load bioinfo-tools java GATK R
		java {params.java_args_gatk} -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T AnalyzeCovariates \
		-R {params.gatk_refpath} \
		-before {input[1]} -after {input[2]} \
		-plots {output[0]}

		java {params.java_args_gatk} -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T PrintReads \
		-R {params.gatk_refpath} \
		-I {input[0]} \
		-BQSR {input[1]} \
		-o {output[1]}
		"""
