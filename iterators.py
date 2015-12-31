#!/usr/bin/env python
import sys
import glob
import pathlib
import json

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

def build_config_dict(sample_name, all_paths):
	config = {
		"samples": sample_name
	}

	config["lanes"] = []
	config["readgroups"] = []	

	for pair in all_paths:
		config["lanes"].append(list(pair))
		config["readgroups"].append(extract_readgroups(pair[0]))
	
	return config

def prettyprint(cfg):
	print(json.dumps(cfg, sort_keys=True, indent=2))

def main():
	paths = all_pairs_for_sample(sys.argv[1])
	cfg = build_config_dict("P1183_101", paths)
	prettyprint(cfg)
	
if __name__ == "__main__":
	main()
