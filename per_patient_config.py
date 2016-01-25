#!/usr/bin/env python
import sys
import json

def main():
	if len(sys.argv) < 2:
		print("{0}: must supple atleast one sample_name.".format(sys.argv[0]))
	else:
		with open("base.json", "r") as handle:
			loaded = json.load(handle)		
		loaded["sample_names"] = sys.argv[1:]
		print(json.dumps(loaded, sort_keys=True, indent=2))

if __name__ == "__main__":
	main()
