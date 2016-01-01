configfile: "config.json"


rule merge_lanes:
	input:
		expand("output/{lane}.bam", lane=config["lanes"].keys())
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
