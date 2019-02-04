import sys

probs = open(sys.argv[1],'r').readlines()
sequence = open(sys.argv[2],'r').readlines()

outf = open(sys.argv[3],'w')

for line in sequence:
    outf.write(line)

for line in probs:
    adist = float(line.split()[0])
    bdist = float(line.split()[1])
    cdist = float(line.split()[2])
    if adist == max(adist,bdist,cdist):
        outf.write('H')
    elif bdist > cdist:
        outf.write('E')
    else:
        outf.write('C')


