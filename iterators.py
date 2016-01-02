#!/usr/bin/env python
import sys
import glob
import pathlib
import os
import json


test_skel = { "reference": "genome.fa" }
raw_data_root = "/proj/b2014316/D.GisselssonNord_15_01/"

def all_pairs_for_sample(path_to_sample):
	to_remove = "1.fastq.gz"
	to_replace = "2.fastq.gz"
	all_paths = []
	for path in glob.glob(path_to_sample+"/*/*"+to_remove):
		p = pathlib.Path(path)
		q = p.with_name(p.name[:-len(to_remove)]+to_replace)
		if q.exists():
			all_paths.append((str(p), str(q)))

	return all_paths

def extract_readgroups(path_to_fq):
	return ""

def build_config_dict(sample_name, all_paths, skeldict={}):
	to_remove = "_1.fastq.gz"
	config = skeldict
	config["samples"] = sample_name
	config["lanes"] = {}

	for pair in all_paths:
		lane_name = os.path.basename(pair[0][:-len(to_remove)])
		config["lanes"][lane_name] = list(pair)
	
	return config

def prettyprint(cfg):
	print(json.dumps(cfg, sort_keys=True, indent=2))

def main():
	if len(sys.argv) < 2:
		print("Usage: {1} [sample_name]", sys.argv[0])
	else:
		path = raw_data_root + sys.argv[1]
		paths = all_pairs_for_sample(path)

		with open("base.json", "r") as handle:
			base_cfg = json.load(handle)	
	
		cfg = build_config_dict(sys.argv[1], paths, skeldict=base_cfg)
		prettyprint(cfg)
	
if __name__ == "__main__":
	main()
