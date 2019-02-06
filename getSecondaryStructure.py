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


def writesspre(filename,outfilename,outfilename2,k):
    infile = open(filename,'r').readlines()
    outfile = open(outfilename,'w')
    o2 = open(outfilename2,'w')
    matcher = readmatcher()
    o2lines = []
    for i,line in enumerate(infile):
        neighbor = ''
        if line.strip() in matcher:
            neighbor = matcher[line.strip()]
            aneighbor = 'XXXXXXXXXXXXXXXXXXX'
            saneighbor = 'CCCCCCCCCCCCCCCCCC'
            adistance = .99
            bneighbor = 'XXXXXXXXXXXXXXXXXXX'
            sbneighbor = 'CCCCCCCCCCCCCCCCCC'
            bdistance = .99
            cneighbor = 'XXXXXXXXXXXXXXXXXXX'
            scneighbor = 'CCCCCCCCCCCCCCCCCC'
            cdistance = .99
            if neighbor[len(neighbor) // 2] == 'A':
                aneighbor = line.strip()
                saneighbor = neighbor
                adistance = .5
            elif neighbor[len(neighbor) // 2] == 'B':
                bneighbor = line.strip()
                sbneighbor = neighbor
                bdistance = .5
            else:
                cneighbor = line.strip()
                scneighbor = neighbor
                cdistance = .5
            for l in range(k):
                outfile.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}\n'.format(i,line[:-1],aneighbor,adistance,bneighbor,bdistance,cneighbor,cdistance,saneighbor,sbneighbor,scneighbor))


        else:
            o2lines.append(line)
    o2.write('# Points: {}\n'.format(len(o2lines) - 1))
    for line in o2lines[1:]:
        o2.write(line)



if len(sys.argv) < 4:    
    writess2(sys.argv[1],sys.argv[2])
else:
    writesspre(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]))
