#!/usr/bin/env python
import sys
import glob
import pathlib
import os
import json
import argparse


raw_data_root = "/proj/b2014316/D.GisselssonNord_15_01/"

# While this assumes that the fastq-files are structured
# according to how the stockholm node of SciLifelab 
def all_pairs_for_sample(path_to_sample):
	to_remove = "1.fastq.gz"
	to_replace = "2.fastq.gz"
	all_paths = []
	# NOTE: recursive globbing using ** and recursive=True requires python >= 3.5
	for path in glob.glob(path_to_sample+"/**/*"+to_remove, recursive=True):
		p = pathlib.Path(path)
		q = p.with_name(p.name[:-len(to_remove)]+to_replace)
		if q.exists():
			all_paths.append((str(p), str(q)))
	return all_paths

def build_config_dict(sample_name, all_paths, config={}):
	to_remove = "_1.fastq.gz"
	
	if "samples" not in config.keys():
		config["samples"] = {}
	
	config["samples"][sample_name] = {}

	for pair in all_paths:
		lane_name = os.path.basename(pair[0][:-len(to_remove)])
		config["samples"][sample_name][lane_name] = list(pair)
	
	return config

def prettyprint(cfg):
	print(json.dumps(cfg, sort_keys=True, indent=2))

def main():
	parser = argparse.ArgumentParser(description='Generate per patient config file for variant calling')
	parser.add_argument("--normal_sample", nargs=1, required=True)
	parser.add_argument("--tumor_samples", nargs='*', required=True)
	parser.add_argument("--dataroot", nargs=1, required=True)
	args = parser.parse_args()
	
	samples = args.normal_sample[:] + args.tumor_samples[:]

	raw_data_root = args.dataroot[0]	

	with open("base.json", "r") as handle:
		cfg = json.load(handle)

	cfg["normal_sample"] = args.normal_sample
	cfg["tumor_sample"] = args.tumor_samples
	cfg["sample_names"] = samples[:]
	
	for sample_path in samples:
		path = raw_data_root + sample_path 
		paths = all_pairs_for_sample(path)

		cfg = build_config_dict(sample_path, paths, config=cfg)
		
	prettyprint(cfg)
	
if __name__ == "__main__":
	main()
