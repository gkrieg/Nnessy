import sys
import traceback
WORDLEN = 23
HALFLEN = WORDLEN // 2
k = 5
numn = 5
k2 = 1000
def getpercentidentity(word,qword):
    pairs = [(word[i],qword[i]) for i in range(len(word))]
    nonzcount = sum([1 if qword[i] != 'Z' else 0 for i in range(len(qword))])
    count = 0
    for l1,l2 in pairs:
        if l1 == l2 and l2 != 'Z':
            count += 1
    return count / nonzcount

def f(x):
    addrange = WORDLEN
    p = .99 #reduction factor
    addfactor = 1
    l = WORDLEN #this is the word length
    z = 2 * (1 - p ** ((addrange+1)/2)) / (1 - p)
    if abs(x) <= addrange // 2:
        if x == 0:
            return (addfactor + 1)/(l+addfactor)
        else:
            return 1/(l+addfactor)
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


def logisticCurve2(distance):
    return 1/distance


def add_probs(probs,i,knn,sstructures,dirnum):
    sstypes = 'GHIEBLTS'
    filename = ''
    bestsstructures = sorted([(knn[a],sstructures[a]) for a in range(len(knn))])
    alphanum = 0
    betanum = 0
    coilnum = 0
    normalizer = 0
    for s in range(len(bestsstructures)):
        #first need to find the value for the normalizer
        if bestsstructures[s][1][HALFLEN] in  'AD':
            alphanum += 1
        elif bestsstructures[s][1][HALFLEN] == 'B':
            betanum += 1
        elif bestsstructures[s][1][HALFLEN] in 'CEFG':
            coilnum += 1
        distancevalue = getdistancevalue(bestsstructures[s][0])
        normalizer += distancevalue

    alphanum = 0
    betanum = 0
    coilnum = 0
    for s in range(len(bestsstructures)):
        if bestsstructures[s][1][HALFLEN] == 'A' or bestsstructures[s][1][HALFLEN] == 'D':
            alphanum += 1
        elif bestsstructures[s][1][HALFLEN] == 'B':
            betanum += 1
        elif bestsstructures[s][1][HALFLEN] in 'CEFG':
            coilnum += 1
        distancevalue = getdistancevalue(bestsstructures[s][0]) / normalizer
        for j in range(i - HALFLEN ,i + HALFLEN + 1):
            relativepos = j - i + HALFLEN
            if j >= 0 and j < len(probs):
                if bestsstructures[s][1][relativepos] in sstypes:
                    sstype = sstypes.index(bestsstructures[s][1][relativepos])
                    probs[j][sstype] += f(j - i) * distancevalue
        

def get_probs(infilename,outfilename,dirnum,ka):
    #this function needs to get all the pieces together from the infile and give them to add_probs
    infile = open(infilename,'r').readlines()
    numn = determinek(infile)
    k = ka
    outfile = open(outfilename,'w')
    identf = open('{}_i'.format(outfilename),'w')
    disf = open('{}_d'.format(outfilename),'w')
    #need to compile the averages here to give to findcoilprobs
    probs = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] for i in range(len(infile) // numn)]
    idents = [[] for i in range(len(infile)//numn)]
    ds = [[] for i in range(len(infile)//numn)]
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
                elif arr[0] == arr[3]:
                    pickedknn = [float(arr[4])]
                elif arr[0] == arr[5]:
                    pickedknn = [float(arr[6])]
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
            add_probs(probs,it//numn,[knns[0]],[sstructures[0]],dirnum)
            peri = getpercentidentity(arr[0].split('_')[0],aas[0])
            dis = neighbors[0][0] 
        for p in range(max(0,(it//numn)-11),min(it//numn + 12,len(infile)//numn)):
            idents[p].append(peri)
            if 'Z' not in infile[it].split()[0]:
                ds[p].append(dis)
    for prob in probs:
        total = sum(prob)
        normprob = [str(p/total) for p in prob]
        outfile.write('{0}\n'.format(' '.join(normprob)))
    for ide in idents:
        identf.write('{}\n'.format(sum(ide)/len(ide)))
    for dis in ds:
        disf.write('{}\n'.format(sum(dis)/len(dis)))
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
                get_probs(inputf,outputf,args.dirnum,args.k)
            except Exception as e:
                print(i)
                traceback.print_exc()

            
process_args()


        
