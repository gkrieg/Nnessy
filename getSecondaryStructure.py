import sys
import pickle

def readmatcher():
    matcherfile = open('matcher.txt','rb')
    matcher = pickle.load(matcherfile)
    return matcher


def writess2(filename,outfilename):
    infile = open(filename,'r').readlines()
    outfile = open(outfilename,'w')
    matcher = readmatcher()
    for line in infile:
        alphaneighbor = line.split()[1]
        betaneighbor = line.split()[3]
        coilneighbor = line.split()[5]
        outfile.write('{0} {1} {2} {3}\n'.format(line[:-1],matcher[alphaneighbor],matcher[betaneighbor],matcher[coilneighbor]))


    
writess2(sys.argv[1],sys.argv[2])
