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

def f(x):
    addrange = WORDLEN
    addfactor = 1
    l = WORDLEN #this is the word length
    if abs(x) <= addrange // 2:
        if x == 0:
            return (addfactor + 1)/(l+addfactor)
        else:
            return 1/(l+addfactor)
    else:
        return 0

def logistic_function(x):
    return 1/(1+np.exp(-x))

def logisticCurve2(distance,filename):
    #read in the logistic curve values
    fi = open(filename,'r').readlines()
    slope = float(fi[0].split()[0][2:-2])
    intercept = float(fi[0].split()[1])
    return logistic_function((distance - intercept) * slope)

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers),1)

def add_probs(probs,i,knn,sstructures,dirnum):
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
    infile = open(infilename,'r').readlines()
    outfile = open(outfilename,'w')
    #need to compile the averages here to give to findcoilprobs
    probs = [[0.0,0.0,0.0] for i in range(len(infile) // numn)]
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
        if sum(alphaknns) <= sum(betaknns) and sum(alphaknns) < sum(coilknns):
            sstructures = alphasstructures
            knns = alphaknns
        elif sum(betaknns) < sum(coilknns):
            sstructures = betasstructures
            knns = betaknns
        else:
            sstructures = coilsstructures
            knns = coilknns
        knns = [float(thing) for thing in knns]
        if ident == False:
            add_probs(probs,it//numn,knns,sstructures,dirnum)
        else:
            add_probs(probs,it//numn,[knns[0]],[sstructures[0]],dirnum)
    for prob in probs:
        total = prob[0] + prob[1] + prob[2]
        outfile.write('{0} {1} {2}\n'.format(prob[0] / total,prob[1] / total,prob[2] / total))
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
    get_probs(args.inputf,args.outputf,args.dirnum,args.k)

            
process_args()


        
