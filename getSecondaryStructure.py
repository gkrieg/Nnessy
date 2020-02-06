import sys
import pickle

def buildmatcher(directorynum):
    alphats = open('/extra/skrieger/code/filegen{}/23alpha.words'.format(directorynum),'r').readlines()
    betats = open('/extra/skrieger/code/filegen{}/23beta.words'.format(directorynum),'r').readlines()
    coilts = open('/extra/skrieger/code/filegen{}/23other.words'.format(directorynum),'r').readlines()

    alphass = open('/extra/skrieger/code/filegen{}/alphaSecondaryStructure.words'.format(directorynum),'r').readlines()
    betass = open('/extra/skrieger/code/filegen{}/betaSecondaryStructure.words'.format(directorynum),'r').readlines()
    coilss = open('/extra/skrieger/code/filegen{}/otherSecondaryStructure.words'.format(directorynum),'r').readlines()

    matcher = {}

    for i in range(len(alphats)):
        matcher[alphats[i].strip()] = alphass[i].strip()
    for i in range(len(betats)):
        matcher[betats[i].strip()] = betass[i].strip()
    for i in range(len(coilts)):
        matcher[coilts[i].strip()] = coilss[i].strip()

    matcherfile = open('matcher{}.txt'.format(directorynum),'wb')
    pickle.dump(matcher,matcherfile)

def piece_together(directorynum):
    alphats = open('../../../code/filegen{}/alpha.words'.format(directorynum),'r').readlines()
    betats = open('../../../code/filegen{}/beta.words'.format(directorynum),'r').readlines()
    coilts = open('../../../code/filegen{}/other.words'.format(directorynum),'r').readlines()

    alphass = open('../../../code/filegen{}/alphasecondaryStructure.txt'.format(directorynum),'r').readlines()
    betass = open('../../../code/filegen{}/betasecondaryStructure.txt'.format(directorynum),'r').readlines()
    coilss = open('../../../code/filegen{}/othersecondaryStructure.txt'.format(directorynum),'r').readlines()

    alphaout = open('../../../code/filegen{}/alphacombined.txt'.format(directorynum),'w')
    betaout = open('../../../code/filegen{}/betacombined.txt'.format(directorynum),'w')
    otherout = open('../../../code/filegen{}/othercombined.txt'.format(directorynum),'w')

    for i in range(len(alphats)):
        alphaout.write('{0} {1}\n'.format(alphats[i].strip(),alphass[i].strip()))
    for i in range(len(betats)):
        betaout.write('{0} {1}\n'.format(betats[i].strip(),betass[i].strip()))
    for i in range(len(coilts)):
        otherout.write('{0} {1}\n'.format(coilts[i].strip(),coilss[i].strip()))


def readmatcher(directorynum):
    matcherfile = open('matcher{}.txt'.format(directorynum),'rb')
    matcher = pickle.load(matcherfile)
    return matcher

def writess(filename):
    infile = open(filename,'r')
    outfile = open('../ss.txt','w')
    matcher = readmatcher()
    for line in infile.readlines():
        outfile.write('{0}\n'.format(matcher[line.split()[1].strip()]))

def writess2(filename,outfilename,directorynum):
    infile = open(filename,'r').readlines()
    outfile = open(outfilename,'w')
    matcher = readmatcher(directorynum)
    for line in infile:
        alphaneighbor = line.split()[1]
        betaneighbor = line.split()[3]
        coilneighbor = line.split()[5]
        outfile.write('{0} {1} {2} {3}\n'.format(line[:-1],matcher[alphaneighbor],matcher[betaneighbor],matcher[coilneighbor]))

def writess3(filename,outfilename):
    infile = open(filename,'r').readlines()
    outfile = open(outfilename,'w')
    matcher = readmatcher()
    for line in infile:
        if len(line.split()) > 1 or len(line.split()) == 0:
            outfile.write(line)
        else:
            outfile.write('{0} {1}\n'.format(line.strip(),matcher[line.strip()]))

def writess4(filename,outfilename):
    infile = open(filename,'r').readlines()
    outfile = open(outfilename,'w')
    matcher = readmatcher()
    for line in infile:
        alphaneighbor = line.split()[1]
        a = line.split()
        betaneighbor = line.split()[4]
        outfile.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(a[0],a[1],a[2],a[4],a[5],matcher[alphaneighbor],matcher[betaneighbor]))
    
def writess5(filename,outfilename,directorynum):
    infile = open(filename,'r').readlines()
    outfile = open(outfilename,'w')
    matcher = readmatcher(directorynum)
    for line in infile:
        s = line.split()
        alphaneighbor = line.split()[1].split('_')[0]
        betaneighbor = line.split()[3].split('_')[0]
        coilneighbor = line.split()[5].split('_')[0]
        newalpha = alphaneighbor + '_' + matcher[alphaneighbor]
        newbeta = betaneighbor + '_' + matcher[betaneighbor]
        newcoil = coilneighbor + '_' + matcher[coilneighbor]
        #qword = s[0].split('_')[0] + '_' + matcher[s[0].split('_')[0]]
        qword = s[0]
        outfile.write('{} {} {} {} {} {} {}\n'.format(qword,newalpha,s[2],newbeta,s[4],newcoil,s[6]))
#buildmatcher(sys.argv[3])
#writess(sys.argv[1])
#writess2(sys.argv[1],sys.argv[2],sys.argv[3])
#writess3(sys.argv[1],sys.argv[2])
#writess4(sys.argv[1],sys.argv[2])
for i in range(int(sys.argv[4])):
    print(i)
    fname = sys.argv[1] + str(i) + '.txt'
    outfname = fname + 'ss'
    try:
        writess5(fname,outfname,sys.argv[3])
    except:
        print('{} did not work'.format(i))
#readmatcher()
#piece_together()
