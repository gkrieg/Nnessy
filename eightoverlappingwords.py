import sys
import traceback
from datetime import datetime
from plotpercentident import getpercentidentity
#from  knn import findcoilprob
#from  knn import findcoilprobs
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
    #p = .8 #reduction factor
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

def determinek(inf):
    first = inf[0].split()[0]
    k = 1
    i = 1
    while inf[i].split()[0] == first:
        k += 1
        i += 1
    return k

def logistic_function(x):
    return 1/(1+np.exp(-x))

def logisticCurve(distance,filename):
    return 1/distance
    ##read in the logistic curve values
    #fi = open(filename,'r').readlines()
    #slope = float(fi[0].split()[0][2:-2])
    #intercept = float(fi[0].split()[1][1:-1])
    #return logistic_function(distance * slope + intercept)

def logisticCurve2(distance,filename):
    return 1/distance
    ##read in the logistic curve values
    #fi = open(filename,'r').readlines()
    #slope = float(fi[0].split()[0][2:-2])
    #intercept = float(fi[0].split()[1])
    #if distance - intercept < 1:
        #return logistic_function((distance - intercept) * slope)
    #else:
        #return .5

def getpoint(y, j,As,notAs):
    margin = 200
    total = 0
    totala = 0
    totalnona = 0
    for i in range(j-margin,j+margin + 1):
        if y[i] == 1:
            #totala += 1.0
            totala += 1.0/As
        else:
            #totalnona += 1.0
            totalnona += 1.0/notAs
        #total += 1 if y[i] == 1 else 0
    return totala / (totala + totalnona)

def meancenter(X,y,As,notAs):
    margin=200
    zerolist = []
    X,y = makeeven(X,y)
    X,y = sortxy(X,y)
    print(X)
    print(y)
    As = len(sum(np.where(y == 1)))
    notAs = len(sum(np.where(y == 0)))
    for i in range(margin,len(X) - margin):
        #print(getpoint(ny,i,As,notAs))
        if abs(getpoint(y,i,As,notAs) - .5) < .01:
            zerolist.append(X[i])
    print(zerolist)
    m = sum(zerolist) / len(zerolist)
    #m = sum(X) / len(X)
    X = X - m
    return X,y,m

def makeeven(X,y):
    zeroes = np.where(y == 0)
    ones = np.where(y == 1)
    smaller = len(zeroes[0]) if len(zeroes[0]) < len(ones[0]) else len(ones[0])
    print('smaller',smaller)
    print('sums',np.sum(y==1))
    print(np.sum(y==0))
    zeroes = np.asarray(random.sample(np.where(y == 0)[0].tolist(),smaller))
    Xtotal = np.append(ones, zeroes)
    X = X[Xtotal]
    y = y[Xtotal]
    

    return X,y

def sortxy(X,y):
    XY = [(X[i],y[i]) for i in range(len(X))]
    XY.sort()
    X = [XY[i][0] for i in range(len(XY))]
    y = [XY[i][1] for i in range(len(XY))]
    return np.asarray(X),np.asarray(y)

def writeLogisticCurves2(filename,outfilename,dirnum):
    inlines = open('{}_{}'.format(filename,dirnum),'r').readlines()
    alphasorted = []
    for line in inlines:
        dist = line.split(',')[0][1:]
        ttype = line.split(',')[1][:-2]
        alphasorted.append((dist,ttype))
    outfile = open('{}.all-{}'.format(outfilename,dirnum),'w')
    alphasorted.sort()
    X = np.asarray([float(x[0]) for x in alphasorted])
    y = np.asarray([int(x[1]) for x in alphasorted])
    X,y,m = meancenter(X,y,0,0)
    X = X.reshape((len(X),1))
    y = y.reshape((len(y),1))
    clf1 = linear_model.LogisticRegression(C=1e5,fit_intercept=False)
    clf1.fit(X,y)
    outfile.write('{0} {1}\n'.format(clf1.coef_,float(m)))




def writeLogisticCurves(filename,outfilename,dirnum):
    inlines = open(filename,'r').readlines()
    alphasorted = []
    for line in inlines:
        adist = line.split(',')[0][1:]
        bdist = line.split(',')[1]
        cdist = line.split(',')[2]
        ttype = line.split(',')[3][:-2]
        alphasorted.append((adist,bdist,cdist,ttype))
    outfile = open('{}.alpha'.format(outfilename),'w')
    outfile2 = open('{}.beta'.format(outfilename),'w')
    outfile3 = open('{}.coil'.format(outfilename),'w')
    X = []
    y1 = []
    y2 = []
    X2 = []
    y3 = []
    X3 = []
    for i,thing in enumerate(alphasorted):
        #if i%1 == 0:
        if float(thing[0]) > .1:
            X.append(float(thing[0]))
            y1.append(1 if thing[3] == 'A' else 0)
        if float(thing[1]) > .1:
            y2.append(1 if thing[3] == 'B' else 0)
            X2.append(float(thing[1]))
        if float(thing[2]) > .1:
            y3.append(1 if thing[3] == 'C' else 0)
            X3.append(float(thing[2]))
    X = np.asarray(X)
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    y3 = np.asarray(y3)
    X3 = np.asarray(X3)
    X2 = np.asarray(X2)
    countA = np.sum(y1 == 1)
    othercountA = len(np.where(y1 == 1)[0])
    print(othercountA,np.sum(y1 == 0))
    countB = np.sum(y2 == 1)
    countC = np.sum(y3 == 1)

    othercountB = len(np.where(y2 == 1)[0])
    print(othercountB,np.sum(y2 == 0))
    As = np.where(y1 == 1)
    Bs = np.where(y2 == 1)
    if countA + countB < countC:
        nonCs = np.where(y3 == 0)
        countnonC = countA + countB
        Cs = np.asarray(random.sample(np.where(y3 == 1)[0].tolist(),countnonC))
    else:
        Cs = np.where(y3 == 1)
        nonCs = np.asarray(random.sample(np.where(y3 == 0)[0].tolist(),countC))
    print('countA and countB',countA,countB)
    nonAs = np.asarray(random.sample(np.where(y1 == 0)[0].tolist(),countA))
    nonBs = np.asarray(random.sample(np.where(y2 == 0)[0].tolist(),countB))
    X,y1,m = meancenter(X,y1,len(nonAs),len(nonAs))
    X2,y2,m2 = meancenter(X2,y2,len(nonBs),len(nonBs))
    X3,y3,m3 = meancenter(X3,y3,len(nonCs),len(nonCs))
    Atotal = np.append(As , nonAs)
    Btotal = np.append(Bs,nonBs)
    Ctotal = np.append(Cs,nonCs)
    print(Atotal)
    #X = X[Atotal]
    #X2 = X2[Btotal]
    #X3 = X3[Ctotal]
    #y1 = y1[Atotal]
    #y2 = y2[Btotal]
    #y3 = y3[Ctotal]
    print('Bs',Bs,len(Bs))
    X = X.reshape((len(X),1))
    X2 = X2.reshape((len(X2),1))
    X3 = X3.reshape((len(X3),1))
    y1 = y1.reshape((len(y1),1))
    y2 = y2.reshape((len(y2),1))
    y3 = y3.reshape((len(y3),1))
    #clf1 = linear_model.LogisticRegression(fit_intercept=True)
    #clf2 = linear_model.LogisticRegression(fit_intercept=True)
    #clf3 = linear_model.LogisticRegression(fit_intercept=True)
    clf1 = linear_model.LogisticRegression(C=1e5,fit_intercept=False)
    clf2 = linear_model.LogisticRegression(C=1e5,fit_intercept=False)
    clf3 = linear_model.LogisticRegression(C=1e5,fit_intercept=False)
    clf1.fit(X,y1)
    clf2.fit(X2,y2)
    clf3.fit(X3,y3)
    print(clf1.coef_,m)
    print(clf2.coef_,m2)
    print(clf3.coef_,m3)
    outfile.write('{0} {1}\n'.format(clf1.coef_,float(m)))
    outfile2.write('{0} {1}\n'.format(clf2.coef_,float(m2)))
    outfile3.write('{0} {1}\n'.format(clf3.coef_,float(m3)))

    


def mean(numbers):
    return float(sum(numbers)) / max(len(numbers),1)

def add_probs(probs,i,knn,sstructures,dirnum):
    sstypes = 'GHIEBLTS'
    filename = ''
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
        if bestsstructures[s][1][HALFLEN] in  'AD':
            filename = 'predictions/{1}/log{0}.alpha'.format(alphanum,dirnum)
            alphanum += 1
        elif bestsstructures[s][1][HALFLEN] == 'B':
            filename = 'predictions/{1}/log{0}.beta'.format(betanum,dirnum)
            betanum += 1
        elif bestsstructures[s][1][HALFLEN] in 'CEFG':
            filename = 'predictions/{1}/log{0}.coil'.format(coilnum,dirnum)
            coilnum += 1
        logisticvalue = logisticCurve2(bestsstructures[s][0],filename)
        normalizer += logisticvalue

    alphanum = 0
    betanum = 0
    coilnum = 0
    for s in range(len(bestsstructures)):
        if bestsstructures[s][1][HALFLEN] == 'A' or bestsstructures[s][1][HALFLEN] == 'D':
            filename = 'predictions/{1}/log{0}.alpha'.format(alphanum,dirnum)
            alphanum += 1
        elif bestsstructures[s][1][HALFLEN] == 'B':
            filename = 'predictions/{1}/log{0}.beta'.format(betanum,dirnum)
            betanum += 1
        elif bestsstructures[s][1][HALFLEN] in 'CEFG':
            filename = 'predictions/{1}/log{0}.coil'.format(coilnum,dirnum)
            coilnum += 1
        logisticvalue = logisticCurve2(bestsstructures[s][0],filename) / normalizer
        for j in range(i - HALFLEN ,i + HALFLEN + 1):
            relativepos = j - i + HALFLEN
            if j >= 0 and j < len(probs):
                if bestsstructures[s][1][relativepos] in sstypes:
                    sstype = sstypes.index(bestsstructures[s][1][relativepos])
                    probs[j][sstype] += f(j - i) * logisticvalue
                #elif bestsstructures[s][1][relativepos] == 'B':
                    #sstype = 1
                    #probs[j][sstype] += f(j - i) * logisticvalue
                #elif bestsstructures[s][1][relativepos] == 'C':
                    #sstype = 2
                    #probs[j][sstype] += f(j - i) * logisticvalue
                #elif bestsstructures[s][1][relativepos] == 'D':
                    #probs[j][0] += .5 * f(j - i) * logisticvalue
                    #probs[j][1] += .5 * f(j - i) * logisticvalue
                #elif bestsstructures[s][1][relativepos] == 'E':
                    #probs[j][0] += .5 * f(j - i) * logisticvalue
                    #probs[j][2] += .5 * f(j - i) * logisticvalue
                #elif bestsstructures[s][1][relativepos] == 'F':
                    #probs[j][1] += .5 * f(j - i) * logisticvalue
                    #probs[j][2] += .5 * f(j - i) * logisticvalue
                #elif bestsstructures[s][1][relativepos] == 'G':
                    #probs[j][0] += .33 * f(j - i) * logisticvalue
                    #probs[j][1] += .33 * f(j - i) * logisticvalue
                    #probs[j][2] += .34 * f(j - i) * logisticvalue
        
def add_probs2(probs,i,knn,sstructures,dirnum):
    #coilfilename = 'predictions/{0}/sortedalpha.txt'.format(dirnum)
    #pc = findcoilprob(mean(knn[:k]),mean(knn[k:]),k2,coilfilename) 
    #probs[i][2] += pc * f(0)
    #then we need to sort the sstructures by the knn distance to just take the closest k of them.
    #bestsstructures = sorted([(knn[a],sstructures[a]) for a in range(len(knn))])[:k]
    bestsstructures = sorted([(knn[a],sstructures[a]) for a in range(len(knn))])
    num = 0
    normalizer = 0
    for s in range(len(bestsstructures)):
        #first need to find the value for the normalizer
        if bestsstructures[s][1][HALFLEN] in  'AD':
            filename = 'predictions/{1}/log{0}.all'.format(num,dirnum)
            num += 1
        elif bestsstructures[s][1][HALFLEN] == 'B':
            filename = 'predictions/{1}/log{0}.all'.format(num,dirnum)
            num += 1
        elif bestsstructures[s][1][HALFLEN] in 'CEFG':
            filename = 'predictions/{1}/log{0}.all'.format(num,dirnum)
            num += 1
        logisticvalue = logisticCurve2(bestsstructures[s][0],filename)
        normalizer += logisticvalue

    num = 0
    for s in range(len(bestsstructures)):
        if bestsstructures[s][1][HALFLEN] == 'A' or bestsstructures[s][1][HALFLEN] == 'D':
            filename = 'predictions/{1}/log{0}.all'.format(num,dirnum)
            num += 1
        elif bestsstructures[s][1][HALFLEN] == 'B':
            filename = 'predictions/{1}/log{0}.all'.format(num,dirnum)
            num += 1
        elif bestsstructures[s][1][HALFLEN] in 'CEFG':
            filename = 'predictions/{1}/log{0}.all'.format(num,dirnum)
            num += 1
        logisticvalue = logisticCurve2(bestsstructures[s][0],filename)
        #logisticvalue = logisticCurve2(bestsstructures[s][0],filename) / normalizer
        for j in range(i - HALFLEN ,i + HALFLEN + 1):
            relativepos = j - i + HALFLEN
            if j >= 0 and j < len(probs):
                if bestsstructures[s][1][relativepos] == 'A':
                    sstype = 0
                    probs[j][sstype] += f(j - i) * logisticvalue
                elif bestsstructures[s][1][relativepos] == 'B':
                    sstype = 1
                    probs[j][sstype] += f(j - i) * logisticvalue
                elif bestsstructures[s][1][relativepos] == 'C':
                    sstype = 2
                    probs[j][sstype] += f(j - i) * logisticvalue
                elif bestsstructures[s][1][relativepos] == 'D':
                    probs[j][0] += .5 * f(j - i) * logisticvalue
                    probs[j][1] += .5 * f(j - i) * logisticvalue
                elif bestsstructures[s][1][relativepos] == 'E':
                    probs[j][0] += .5 * f(j - i) * logisticvalue
                    probs[j][2] += .5 * f(j - i) * logisticvalue
                elif bestsstructures[s][1][relativepos] == 'F':
                    probs[j][1] += .5 * f(j - i) * logisticvalue
                    probs[j][2] += .5 * f(j - i) * logisticvalue
                elif bestsstructures[s][1][relativepos] == 'G':
                    probs[j][0] += .33 * f(j - i) * logisticvalue
                    probs[j][1] += .33 * f(j - i) * logisticvalue
                    probs[j][2] += .34 * f(j - i) * logisticvalue

def get_probs(infilename,outfilename,dirnum,ka):
    #this function needs to get all the pieces together from the infile and give them to add_probs
    print(infilename)
    print('before reading file',str(datetime.now().time()))
    infile = open(infilename,'r').readlines()
    numn = determinek(infile)
    k = ka
    #print('after reading file',str(datetime.now().time()))
    outfile = open(outfilename,'w')
    identf = open('{}_i'.format(outfilename),'w')
    disf = open('{}_d'.format(outfilename),'w')
    #need to compile the averages here to give to findcoilprobs
    probs = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] for i in range(len(infile) // numn)]
    idents = [[] for i in range(len(infile)//numn)]
    ds = [[] for i in range(len(infile)//numn)]
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
        pickedknn = []
        pickedss = []
        alphaaa = []
        betaaa = []
        coilaa = []
        for num in range(it,it+k):
            arr = infile[num].split()
            if arr[0] == arr[1] or arr[0] == arr[3] or arr[0] == arr[5]:
                ident = True
                if arr[0] == arr[1]:
                    pickedknn = [float(arr[2])]
                    #pickedss = [arr[7]]
                elif arr[0] == arr[3]:
                    pickedknn = [float(arr[4])]
                    #pickedss = [arr[8]]
                elif arr[0] == arr[5]:
                    pickedknn = [float(arr[6])]
                    #pickedss = [arr[9]]
            alphaknns.append(float(arr[2]))
            betaknns.append(float(arr[4]))
            coilknns.append(float(arr[6]))
            alphasstructures.append(arr[1].split('_')[1])
            betasstructures.append(arr[3].split('_')[1])
            coilsstructures.append(arr[5].split('_')[1])
            alphaaa.append(arr[1].split('_')[0])
            betaaa.append(arr[3].split('_')[0])
            coilaa.append(arr[5].split('_')[0])

        neighbors = sorted([(alphaknns[i],alphasstructures[i],alphaaa[i]) for i in range(len(alphaknns))] + [(betaknns[j],betasstructures[j],betaaa[j]) for j in range(len(betaknns))] + [(coilknns[l],coilsstructures[l],coilaa[l]) for l in range(len(coilknns))])
        #This section for removing identities only
        #ident == False
        #for e,n in enumerate(neighbors):
            #if n[2] == infile[it].split()[0].split('_')[0]:
                #del neighbors[e]
        #end section
        knns = [neighbors[i][0] for i in range(k)]
        sstructures = [neighbors[i][1] for i in range(k)]
        aas = [neighbors[i][2] for i in range(k)]
        if ident == False:
            add_probs(probs,it//numn,knns,sstructures,dirnum)
            peri = 0
            for aa in aas:
                peri += getpercentidentity(arr[0].split('_')[0],aa)
            peri /= len(aas)
            dis = sum([neighbors[q][0] for q in range(k)])/k
        else:
            #add_probs(probs,it//numn,pickedknn,pickedss,dirnum)
            add_probs(probs,it//numn,[knns[0]],[sstructures[0]],dirnum)
            peri = getpercentidentity(arr[0].split('_')[0],aas[0])
            dis = neighbors[0][0] 
        for p in range(max(0,(it//numn)-11),min(it//numn + 12,len(infile)//numn)):
            idents[p].append(peri)
            if 'Z' not in infile[it].split()[0]:
                ds[p].append(dis)
    print('after adding probs',str(datetime.now().time()))
    for prob in probs:
        #total = prob[0] + prob[1] + prob[2]
        #outfile.write('{0} {1} {2}\n'.format(prob[0] / total,prob[1] / total,prob[2] / total))
        #use the following two lines for the different kind of normalization
        total = sum(prob)
        normprob = [str(p/total) for p in prob]
        outfile.write('{0}\n'.format(' '.join(normprob)))
    for ide in idents:
        identf.write('{}\n'.format(sum(ide)/len(ide)))
    for dis in ds:
        disf.write('{}\n'.format(sum(dis)/len(dis)))
    print('after printing probs',str(datetime.now().time()))
    return probs

def process_args():
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Parses the input.')
    parser.add_argument("inputf", help="the input filename for which you want the class-membership probabilities", type=str, default="predictions/12222/predictknn1.txtss")
    parser.add_argument("outputf", help="the output filename to which you want to write the class-membership probabilities", type=str, default="outtt")
    parser.add_argument("dirnum", help="the directory number", type=str, default='1')
    parser.add_argument("k", help="the number of neighbors to look at", type=int,default=3)
    parser.add_argument("count", type=int, default=1)
    parser.add_argument("-l", help="flag for generating logistic curves", action='store_true')
    args = parser.parse_args()
    print('k = {}'.format(args.k))
    for i in range(args.count):
        inputf = '{}{}.txtss'.format(args.inputf,i)
        outputf = '{}{}'.format(args.outputf,i)
        if args.l == True:
            writeLogisticCurves2(args.inputf,args.outputf,args.dirnum)
        else:
            try:
                #get_probs_targets(inputf,outputf,args.dirnum,args.k)
                get_probs(inputf,outputf,args.dirnum,args.k)
            #get_probs_ident(inputf,outputf,args.dirnum,args.k)
            except Exception as e:
                print(i)
                traceback.print_exc()

            
#get_probs('predictions/12222/predictknn1.txtss','outtt')
#writeLogisticCurves('coildistances2/sortedalpha.txt')
process_args()


        
