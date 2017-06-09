#!/usr/bin/python

f = open("5GEM_mapping.cfg", "r")

lastplane = "u0"

for line in f:
    if line[0] != '#':
	vals = line.replace('\n', '').replace('\t', '').split(',')
	if len(vals) > 1:
	  orient = 'x'
	  if vals[2] == '0':
	      orient = 'x'
	  else:
	      orient = 'y'
	  gemname =  orient + str(int(vals[1])+1)
	  if gemname != lastplane:
 	      print "sbs.gems." + gemname + ".chanmap = \\"
	      lastplane = gemname
	  chanline = "11 0 " + vals[0] + " " + vals[1] + " " + \
	      vals[3] + " " + vals[4] + " " + vals[5] + " " + \
	      vals[6] + " \\"
	  print chanline
	

