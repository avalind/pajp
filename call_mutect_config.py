#!/usr/bin/env python
import sys
import json
import argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--normal_sample", nargs=1, required=True)
	parser.add_argument("--tumor_samples", nargs='*', required=True)
	args = parser.parse_args()

	with open("base.json", "r") as handle:
		loaded = json.load(handle)

	loaded["normal_sample"] = args.normal_sample
	loaded["tumor_sample"] = args.tumor_samples
	loaded["sample_names"] = args.tumor_samples[:] + args.normal_sample[:]
	print(json.dumps(loaded, sort_keys=True, indent=2))

if __name__ == "__main__":
	main()
