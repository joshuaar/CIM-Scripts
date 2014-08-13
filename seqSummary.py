#!/usr/bin/python
import argparse,Data,csv,sys
import numpy as np

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
        parser.add_argument("infile")
        args = parser.parse_args()

	x = Data.ChipData.fromFileName(args.infile)
	
	mins = []
	maxs = []
	means = []
	rngs = []
	
	for i in x:
		vals = i.values.transpose()[0]
		mean = np.mean(vals)
		mx = max(vals)
		mn = min(vals)
		rng = mx - mn	
		mins.append(mn)
		maxs.append(mx)
		means.append(mean)
		rngs.append(rng)
	totals = zip(x.samples,mins,maxs,rngs,means)	
	totals = [["Sample","Min","Max","Range","Mean"]]+totals
	wtr = csv.writer(sys.stdout)
	wtr.writerows(totals)
	
	
	
	
