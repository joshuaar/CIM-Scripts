#!/usr/bin/python
import sys
if __name__ == "__main__":
	infile = sys.argv[1]
	lst = open(infile).read().split("\n")
	for i in lst:
		if len(i) > 0:
			print ">{0}".format(i)
			print i

