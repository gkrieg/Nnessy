import sys
import pickle


def readmatcher(directorynum):
    matcherfile = open('matcher{}.txt'.format(directorynum),'rb')
    matcher = pickle.load(matcherfile)
    return matcher

    
def writess8(filename,outfilename,directorynum):
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
        qword = s[0]
        outfile.write('{} {} {} {} {} {} {}\n'.format(qword,newalpha,s[2],newbeta,s[4],newcoil,s[6]))

for i in range(int(sys.argv[4])):
    print(i)
    fname = sys.argv[1]
    outfname = sys.argv[2]
    try:
        writess8(fname,outfname,sys.argv[3])
    except:
        print('{} did not work'.format(fname))
