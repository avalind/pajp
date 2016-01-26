#!/usr/bin/env python
import sys
import json
import argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("normal_sample")
	parser.add_argument("tumor_sample")
	args = parser.parse_args()

	with open("base.json", "r") as handle:
		loaded = json.load(handle)
	
	loaded["normal_sample"] = args.normal_sample
	loaded["tumor_sample"] = args.tumor_sample
	print(json.dumps(loaded, sort_keys=True, indent=2))

if __name__ == "__main__":
	main()
