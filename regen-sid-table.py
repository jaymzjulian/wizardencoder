#!/bim/python
FREQSCALE=3

v1=[]
v2=[]

for val in range(256):
         v1.append(((val)<<(8-(FREQSCALE-1)))&0xff)
         v2.append(((val)>>((FREQSCALE-1)))&0xff)
print "sflo"
for val in range(256):
	print ".byte",v1[val]
print "sfhi"
for val in range(256):
	print ".byte",v2[val]
