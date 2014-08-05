#!/usr/bin/python
import gzip,csv,sys,argparse,os
"Turns a gpr file into a row of data. Works."

def readMeta(handle):
	"reads GPR metadata, returns a csv reader object and the metadata as a dict"
	meta = {}
	rdr = csv.reader(handle)
	version = rdr.next()[0]
	headers,cols = [int(i) for i in rdr.next()[0].replace(" ","").split("\t")]
	meta = {j[0]:j[1] for j in [rdr.next()[0].split("=") for i in range(headers)]}
	return rdr,meta

def readGPR(rdr,col):
	"Extracts a column, matched to the col argument, and prints to stdout as a row"
	colNames = rdr.next()[0].replace('"',"").replace("'","").split("\t")
	colOfInterest = colNames.index(col)
	sys.stdout.write(name)
	for i in rdr:
		sys.stdout.write("\t"+i[0].replace('"',"").split("\t")[colOfInterest])
	sys.stdout.write("\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Extract a column of GPR file and print to stdout')
	parser.add_argument("--gzip", action="store_true")	
	parser.add_argument("--jpeg", action="store_true")	
	parser.add_argument("--col", nargs=1)
	parser.add_argument("input", nargs=1)
	args = parser.parse_args()
	if args.gzip:
		fHandle = gzip.open(args.input[0])
	else:
		fHandle = open(args.input)
	rdr,meta = readMeta(fHandle)
	if args.jpeg:
		name = os.path.basename(meta["JpegImage"])
	else:
		name = args.input
	sys.stdout.write(name)
	readGPR(rdr,args.col[0])
	
	
	
	

		
