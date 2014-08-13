import sys
if __name__ == "__main":
	infile = sys.argv[1]
	lst = open(infile).read().split("\n")
	for i in lst:
		print ">{0}".format(i)
		print i

