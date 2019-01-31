import sys
from datetime import datetime
#from sklearn import linear_model
#from scipy import optimize
import random
#import numpy as np
WORDLEN = 23
HALFLEN = WORDLEN // 2
k = 5
numn = 5
k2 = 1000

#def f(x):
    #p = .3 #reduction factor
    ##p = .7 #reduction factor
    #l = 23 #this is the word length
    #z = 2 * (1 - p ** ((l+1)/2)) / (1 - p)
    #if x == 0:
        #return 2 * p**(abs(x)) / z
    #else:
        #return p**(abs(x)) / z
def f(x):
    addrange = WORDLEN
    a  = [0.8230912476722533, 0.827735774657546, 0.8315047355992601, 0.8363320184925987, 0.8406575314708341, 0.8449112013246257, 0.8489696198930714, 0.8531684183329027, 0.8586477909044082, 0.86688338613106, 0.8769002626194007, 0.8924968108546909, 0.8803100954968237, 0.8699236843631664, 0.8616165890816758, 0.8558733800510403, 0.851171279281982, 0.8477056400056826, 0.8431070174409346, 0.8390498213159554, 0.8335189591454463, 0.8271263083020018, 0.8204840573798354]
    #p = 0.993472336883 #reduction factor
    p = .99 #reduction factor
    addfactor = 1
    #p = .7 #reduction factor
    l = WORDLEN #this is the word length
    z = 2 * (1 - p ** ((addrange+1)/2)) / (1 - p)
    if abs(x) <= addrange // 2:
        if x == 0:
            #return 2 * p**(abs(x)) / z
            return (addfactor + 1)/(l+addfactor)
            #return 1
        else:
            #return p**(abs(x)) / z
            return 1/(l+addfactor)
            #return 0
    else:
        return 0

def logistic_function(x):
    return 1/(1+np.exp(-x))

def logisticCurve(distance,filename):
    #read in the logistic curve values
    fi = open(filename,'r').readlines()
    slope = float(fi[0].split()[0][2:-2])
    intercept = float(fi[0].split()[1][1:-1])
    return logistic_function(distance * slope + intercept)

def logisticCurve2(distance,filename):
    #read in the logistic curve values
    fi = open(filename,'r').readlines()
    slope = float(fi[0].split()[0][2:-2])
    intercept = float(fi[0].split()[1])
    return logistic_function((distance - intercept) * slope)

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers),1)

def add_probs(probs,i,knn,sstructures,dirnum):
    #coilfilename = 'predictions/{0}/sortedalpha.txt'.format(dirnum)
    #pc = findcoilprob(mean(knn[:k]),mean(knn[k:]),k2,coilfilename) 
    #probs[i][2] += pc * f(0)
    #then we need to sort the sstructures by the knn distance to just take the closest k of them.
    #bestsstructures = sorted([(knn[a],sstructures[a]) for a in range(len(knn))])[:k]
    bestsstructures = sorted([(knn[a],sstructures[a]) for a in range(len(knn))])
    alphanum = 0
    betanum = 0
    coilnum = 0
    normalizer = 0
    for s in range(len(bestsstructures)):
        #first need to find the value for the normalizer
        if bestsstructures[s][1][HALFLEN] == 'A':
            filename = 'logs/log{0}.alpha'.format(alphanum)
            alphanum += 1
        elif bestsstructures[s][1][HALFLEN] == 'B':
            filename = 'logs/log{0}.beta'.format(betanum)
            betanum += 1
        elif bestsstructures[s][1][HALFLEN] == 'C':
            filename = 'logs/log{0}.coil'.format(coilnum)
            coilnum += 1
        logisticvalue = logisticCurve2(bestsstructures[s][0],filename)
        normalizer += logisticvalue

    alphanum = 0
    betanum = 0
    coilnum = 0
    for s in range(len(bestsstructures)):
        if bestsstructures[s][1][HALFLEN] == 'A':
            filename = 'logs/log{0}.alpha'.format(alphanum)
            alphanum += 1
        elif bestsstructures[s][1][HALFLEN] == 'B':
            filename = 'logs/log{0}.beta'.format(betanum)
            betanum += 1
        elif bestsstructures[s][1][HALFLEN] == 'C':
            filename = 'logs/log{0}.coil'.format(coilnum)
            coilnum += 1
        logisticvalue = logisticCurve2(bestsstructures[s][0],filename) / normalizer
        for j in range(i - HALFLEN ,i + HALFLEN + 1):
            relativepos = j - i + HALFLEN
            if j >= 0 and j < len(probs):
                if bestsstructures[s][1][relativepos] == 'A':
                    sstype = 0
                elif bestsstructures[s][1][relativepos] == 'B':
                    sstype = 1
                elif bestsstructures[s][1][relativepos] == 'C':
                    sstype = 2
                probs[j][sstype] += f(j - i) * logisticvalue
        

def get_probs(infilename,outfilename,dirnum,ka):
    #this function needs to get all the pieces together from the infile and give them to add_probs
    k = ka
    print(infilename)
    print('before reading file',str(datetime.now().time()))
    infile = open(infilename,'r').readlines()
    #print('after reading file',str(datetime.now().time()))
    outfile = open(outfilename,'w')
    #need to compile the averages here to give to findcoilprobs
    probs = [[0.0,0.0,0.0] for i in range(len(infile) // numn)]
    print(len(probs),'len probs')
    #alphadists = []
    #betadists = []
    #coildists = []
    #for it in range(0,len(infile),numn):
        #alphaknns = []
        #betaknns = []
        #coilknns = []
        #for num in range(it,min(it+k,len(infile))):
            #arr = infile[num].split()
            #alphaknns.append(arr[2])
            #betaknns.append(arr[4])
            #coilknns.append(arr[6])
        #knns = alphaknns + betaknns + coilknns
        #knns = [float(thing) for thing in knns]
        #alphadists.append(mean(knns[:k]))
        #betadists.append(mean(knns[k:]))
        #coildists.append(mean(knns[k:]))
    #print('after getting lists',str(datetime.now().time()))
    #coilfilename = 'predictions/{1}/sortedalpha.txt{0}'.format(k,dirnum)
    #coilfilename = 'predictions/{1}/sortedalpha.txt{0}'.format(0,dirnum)
    #print(coilfilename)
    #coilprobs = findcoilprobs(alphadists,betadists,k2,coilfilename)
    #print('after getting coilprobs',str(datetime.now().time()))
    for it in range(0,len(infile),numn):
        alphaknns = []
        betaknns = []
        coilknns = []
        alphasstructures = []
        betasstructures = []
        coilsstructures = []
        ident = False
        for num in range(it,it+k):
            arr = infile[num].split()
            if arr[0] == arr[1] or arr[0] == arr[3] or arr[0] == arr[5]:
                ident = True
            alphaknns.append(float(arr[2]))
            betaknns.append(float(arr[4]))
            coilknns.append(float(arr[6]))
            alphasstructures.append(arr[7])
            betasstructures.append(arr[8])
            coilsstructures.append(arr[9])
        #knns = alphaknns + betaknns + coilknns
        if sum(alphaknns) <= sum(betaknns) and sum(alphaknns) < sum(coilknns):
            #sstructures = [alphasstructures[0]]
            #knns = [alphaknns[0]]
            sstructures = alphasstructures
            knns = alphaknns
        elif sum(betaknns) < sum(coilknns):
            #sstructures = [betasstructures[0]]
            #knns = [betaknns[0]]
            sstructures = betasstructures
            knns = betaknns
        else:
            #sstructures = [coilsstructures[0]]
            #knns = [coilknns[0]]
            sstructures = coilsstructures
            knns = coilknns
        knns = [float(thing) for thing in knns]
        #sstructures = alphasstructures + betasstructures + coilsstructures
        if ident == False:
            add_probs(probs,it//numn,knns,sstructures,dirnum)
        else:
            add_probs(probs,it//numn,[knns[0]],[sstructures[0]],dirnum)
    print('after adding probs',str(datetime.now().time()))
    for prob in probs:
        #total = prob[0] + prob[1] + prob[2]
        #outfile.write('{0} {1} {2}\n'.format(prob[0] / total,prob[1] / total,prob[2] / total))
        #use the following two lines for the different kind of normalization
        total = prob[0] + prob[1] + prob[2]
        outfile.write('{0} {1} {2}\n'.format(prob[0] / total,prob[1] / total,prob[2] / total))
    print('after printing probs',str(datetime.now().time()))
    return probs

def process_args():
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Parses the input.')
    parser.add_argument("inputf", help="the input filename for which you want the class-membership probabilities", type=str, default="predictions/12222/predictknn1.txtss")
    parser.add_argument("outputf", help="the output filename to which you want to write the class-membership probabilities", type=str, default="outtt")
    parser.add_argument("dirnum", help="the directory number", type=str, default='1')
    parser.add_argument("k", help="the number of neighbors to look at", type=int,default=5)
    args = parser.parse_args()
    print('k = {}'.format(args.k))
    get_probs(args.inputf,args.outputf,args.dirnum,args.k)

            
#get_probs('predictions/12222/predictknn1.txtss','outtt')
#writeLogisticCurves('coildistances2/sortedalpha.txt')
process_args()


        
