#!/bin/python
import math

cfile=open("sidc.h", "wt")

root=pow(2, 1.0/36.0)
# thanks codebase!
sidconstant=17.02841924

base=440.0/32.0
print "Base Freq:",base
print "Root",root

print >>cfile, "static int sidFreq[256]={"
for c in range(256):
	final=base*pow(root, c)
	print c,final,final*sidconstant
	print >>cfile, int(final*sidconstant),","

print >>cfile, "0};"
