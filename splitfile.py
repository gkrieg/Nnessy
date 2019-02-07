from sys import argv

infile = open(argv[1],'r').readlines()
numsplits = int(argv[2])

fileit = 1
splitsize = (len(infile) - 1) // numsplits
leftover = (len(infile) - 1 ) % numsplits
i = 1
n = 0
print('len',len(infile))
while i < len(infile):
    print(i)
    outfile = open(argv[1] + "_" + str(numsplits) + "_" + str(n),'w')
    numpoints = splitsize
    if n < leftover:
        numpoints += 1
    outfile.write('# Points: {}\n'.format(numpoints))
    for j in range(i,i + numpoints):
        outfile.write(infile[j])
    outfile.close()
    i += numpoints
    n += 1


