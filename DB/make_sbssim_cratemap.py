#!/usr/bin/env python
import re ## Regular expressions

known_adcs = [ '250',  '792', '1881' ]
known_tdcs = ['1190', '1877', '3204', '1875', '3201', '775']
known_mpds = ['3561']
known_scalers = [ '1151', '3800', '3801', '560' ]

def getSimEncoder(val):
    for adc in known_adcs:
        if val == adc:
            return '50250'
    for tdc in known_tdcs:
        if val == tdc:
            return '53240'
    for mpd in known_mpds:
        if val == mpd:
            return '53561'
    return val ## encoder not found, good luck!

def makeSimMap():
    f = open("db_cratemap.dat",'r');
    for line in f:
        line = line.rstrip() ## Remove newline character
        line_split = line.split('#',1) ## Split between comments and not
        ## Just print empty lines, or lines with only comments
        ## all without modification
        if len(line_split) == 0 or not line_split[0].strip():
            print(line)
        ## Print the cratemap lines
        elif re.match(r'^====',line_split[0].strip()):
            print(line)
        else:
            values = line_split[0].strip().split() ## Split on space
            if len(values) > 1: ## Good, at least the module number is specified
                values[1] = getSimEncoder(values[1])

            newline = ''
            for val in values:
                newline = newline+'%6s '%val
            if len(line_split) > 1:
                newline = newline+' #'+line_split[1]
            print(newline)

        #print('%d :: %s'%(len(line_split),line))
        #elif len(line_split) 1: ## Only a comment
        #elif re.match(r'^#',line):
#            print(line)
        #items = line.split('#',1) ## Split out comments
        ## Preserve comments as is
#        if re.match(r'^====',line): ## Preserve crate map definition
#            print(line)
#        else: ## Anything else should be the actual crate map
#            items = line.split()



makeSimMap()
