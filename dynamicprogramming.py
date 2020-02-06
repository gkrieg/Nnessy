
from numpy import inf
from math import log

def getprxati(i,rand):
    prxati = []
    neginf = -1000000
    fin = open(rand,'r').readlines()
    #print(len(fin))
    #fin = open('../lpgeneration/predictions/' + rand + '/newprobs' + str(i),'r').readlines()
    #fin = open('../lpgeneration/outtt','r').readlines()
    for line in fin:
        single = []
        needschange = False
        for i,word in enumerate(line.strip().split()):
            if word < '0.00001':
                single.append(neginf)
            else:
                single.append(log(float(word)))
            if float(word) > .99:
                needschange = True
        if needschange == True:
            single2 = []
            for i,thing in enumerate(single):
                if thing > log(.99):
                    single2.append(log(1-.01*len(single)))
                else:
                    single2.append(log(.01))
            prxati.append(single2)

        else:
            prxati.append(single)
    return prxati
        
def highestprobability(prxati):
    answer = []
    for line in prxati:
        if line[0] >= line[1] and line[0] >= line[2]:
            answer.append('A')
        elif line[1] >= line[2]:
            answer.append('B')
        else:
            answer.append('C')
    return ''.join(answer)


def classfind(a):
    if a == 'G' or a == 'H' or a == 'I' or a == 'A':
        return 'A'
    elif a == 'B' or a == 'E' or a == 'B':
        return 'B'
    else:
        return 'C'
def gettransitions():
    #lines = open('../../code/filegen4/sequences/newpdb.seqs','r').readlines()
    lines = open('cleanpdb.seqs','r').readlines()
    maxlength = 200
    counts = [[0 for j in range(maxlength)] for i in range(9)]
    transitions = [[0 for j in range(maxlength)] for i in range(9)]
    #for i in range(3,len(lines),5):
    for i in range(1,len(lines),2):
        line = lines[i].strip()
        length = 0
        for j in range(len(line)-1):
            if length > maxlength - 1:
                length = maxlength - 1
            currchar = classfind(line[j])
            nextchar = classfind(line[j+1])
            if currchar == 'A':
                if nextchar == 'A':
                    counts[0][length] += 1
                    length += 1
                if nextchar == 'B':
                    counts[1][length] += 1
                    length = 0
                if nextchar == 'C':
                    counts[2][length] += 1
                    length = 0
            if currchar == 'B':
                if nextchar == 'A':
                    counts[3][length] += 1
                    length = 0
                if nextchar == 'B':
                    counts[4][length] += 1
                    length += 1
                if nextchar == 'C':
                    counts[5][length] += 1
                    length = 0
            if currchar == 'C':
                if nextchar == 'A':
                    counts[6][length] += 1
                    length = 0
                if nextchar == 'B':
                    counts[7][length] += 1
                    length = 0
                if nextchar == 'C':
                    counts[8][length] += 1
                    length += 1
    #print(len(counts[0]),len(counts))
    #print(counts)
    #print(len(transitions[0]),len(transitions))
    for i in range(len(counts[0])):
        acounts = sum([counts[j][i] for j in range(0,3)])
        bcounts = sum([counts[j][i] for j in range(3,6)])
        ccounts = sum([counts[j][i] for j in range(6,9)])
        transitions[0][i] = log(counts[0][i] / acounts) if counts[0][i] != 0 and acounts != 0 else -10000000
        transitions[1][i] = log(counts[1][i] / acounts) if counts[1][i] != 0  and acounts != 0 else -10000000
        transitions[2][i] = log(counts[2][i] / acounts) if counts[2][i] != 0 and acounts != 0 else -10000000
        transitions[3][i] = log(counts[3][i] / bcounts) if counts[3][i] != 0 and bcounts != 0 else -10000000
        transitions[4][i] = log(counts[4][i] / bcounts) if counts[4][i] != 0 and bcounts != 0 else -10000000
        transitions[5][i] = log(counts[5][i] / bcounts) if counts[5][i] != 0 and bcounts != 0 else -10000000
        transitions[6][i] = log(counts[6][i] / ccounts) if counts[6][i] != 0 and ccounts != 0 else -10000000
        transitions[7][i] = log(counts[7][i] / ccounts) if counts[7][i] != 0 and ccounts != 0 else -10000000
        transitions[8][i] = log(counts[8][i] / ccounts) if counts[8][i] != 0 and ccounts != 0 else -10000000
    return transitions
 
def gettransitions8():
    lines = open('/extra/skrieger/code/filegen4/sequences/TIMELINE/newpdb.seqs','r').readlines()
    #lines = open('cleanpdb.seqs','r').readlines()
    maxlength = 200
    letters = ['G','H','I','E','B','L','T','S']
    counts = [[[0 for j in range(maxlength)] for i in range(8)] for k in range(8)]
    transitions = [[[0 for j in range(maxlength)] for i in range(8)] for k in range(8)]
    for i in range(3,len(lines),5):
    #for i in range(1,len(lines),2):
        line = lines[i].strip()
        length = 0
        for j in range(len(line)-1):
            if length > maxlength - 1:
                length = maxlength - 1
            currchar = line[j]
            nextchar = line[j+1]
            currnum = letters.index(currchar)
            nextnum = letters.index(nextchar)
            if currnum == nextnum:
                counts[currnum][nextnum][length] += 1
                length += 1
            else:
                counts[currnum][nextnum][length] += 1
                counts[nextnum][currnum][length] += 1
                length = 0
    for i in range(len(counts)):
        for j in range(len(counts[i][i])):
            acounts = (sum([counts[i][k][j] for k in range(len(counts))]) + counts[i][i][j]) // 2
            for k in range(len(counts[i])):
                transitions[i][k][j] = log(counts[i][k][j] / acounts) if counts[i][k][j] != 0 and acounts != 0 else -10000000
    #transitions8 = [[[0 for a in range(maxlength)] for b in range(8)] for c in range(8)]
    #for i in range(maxlength):
        #for k in range(8):
            #if k < 3:
                #transitions8[k][0][i] = transitions[0][i]
                #transitions8[k][1][i] = transitions[0][i]
                #transitions8[k][2][i] = transitions[0][i]
                #transitions8[0][k][i] = transitions[0][i]
                #transitions8[1][k][i] = transitions[0][i]
                #transitions8[2][k][i] = transitions[0][i]
#
                #transitions8[k][3][i] = transitions[1][i]
                #transitions8[k][4][i] = transitions[1][i]
                #transitions8[3][k][i] = transitions[1][i]
                #transitions8[4][k][i] = transitions[1][i]
#
                #transitions8[k][5][i] = transitions[2][i]
                #transitions8[k][6][i] = transitions[2][i]
                #transitions8[k][7][i] = transitions[2][i]
                #transitions8[5][k][i] = transitions[2][i]
                #transitions8[6][k][i] = transitions[2][i]
                #transitions8[7][k][i] = transitions[2][i]

            #elif k < 5:
                #transitions8[k][0][i] = transitions[3][i]
                #transitions8[k][1][i] = transitions[3][i]
                #transitions8[k][2][i] = transitions[3][i]
                #transitions8[0][k][i] = transitions[3][i]
                #transitions8[1][k][i] = transitions[3][i]
                #transitions8[2][k][i] = transitions[3][i]
#
                #transitions8[k][3][i] = transitions[4][i]
                #transitions8[k][4][i] = transitions[4][i]
                #transitions8[3][k][i] = transitions[4][i]
                #transitions8[4][k][i] = transitions[4][i]
#
                #transitions8[k][5][i] = transitions[5][i]
                #transitions8[k][6][i] = transitions[5][i]
                #transitions8[k][7][i] = transitions[5][i]
                #transitions8[5][k][i] = transitions[5][i]
                #transitions8[6][k][i] = transitions[5][i]
                #transitions8[7][k][i] = transitions[5][i]
            #else:
                #transitions8[k][0][i] = transitions[6][i]
                #transitions8[k][1][i] = transitions[6][i]
                #transitions8[k][2][i] = transitions[6][i]
                #transitions8[0][k][i] = transitions[6][i]
                #transitions8[1][k][i] = transitions[6][i]
                #transitions8[2][k][i] = transitions[6][i]
#
                #transitions8[k][3][i] = transitions[7][i]
                #transitions8[k][4][i] = transitions[7][i]
                #transitions8[3][k][i] = transitions[7][i]
                #transitions8[4][k][i] = transitions[7][i]
#
                #transitions8[k][5][i] = transitions[8][i]
                #transitions8[k][6][i] = transitions[8][i]
                #transitions8[k][7][i] = transitions[8][i]
                #transitions8[5][k][i] = transitions[8][i]
                #transitions8[6][k][i] = transitions[8][i]
                #transitions8[7][k][i] = transitions[8][i]
#
    return transitions


def gettransitions2():

    maxlength = 500
    transitions = [[0 for j in range(maxlength)] for i in range(9)]
    for i in range(len(transitions[0])):
        transitions[0][i] = 0
        transitions[1][i] = log(.01)
        transitions[2][i] =  log(.01)
        transitions[3][i] =  log(.01)
        transitions[4][i] = 0
        transitions[5][i] =  log(.01)
        transitions[6][i] =  log(.01)
        transitions[7][i] =  log(.01)
        transitions[8][i] = 0
    return transitions
        


def fillmatrix():
    global prxati
    transcoefficient = 1
    maxAhistory = 80
    maxBhistory = 30
    maxChistory = 40
    minArun = 3
    minBrun = 1
    minCrun = 1
    import sys
    w = float(sys.argv[4])
    x = float(sys.argv[2])
    y = float(sys.argv[3])
    y = 0
    z = 1
    prarun = [1e-50, 1e-50, 0.1802433115377602, 0.09596130634503, 0.04924410173583071, 0.0528738027630882, 0.05563515426229366, 0.05408949210236104, 0.056451402818662566, 0.0616788669775358, 0.05791023002578998, 0.04990404737715025, 0.05156259497572963, 0.044581064770191296, 0.03272809371228107, 0.02647597711030644, 0.02375804308750358, 0.019650749819817473, 0.01354625263761169, 0.011800870086227109, 0.010046804039562004, 0.008188536049530657, 0.005913460286034335, 0.005679005913460286, 0.004480683564748482, 0.003673118503660093, 0.0036991689895016543, 0.0037425864659042557, 0.0026571495558392164, 0.002153506829569038, 0.0019972039145196726, 0.0012243728345533644, 0.0007207301082831862, 0.0010941204053455597, 0.0007120466130026658, 0.0007815145752468283, 0.0007467805941247471, 0.0005470602026727799, 0.000494959230989658]
    #prarun = [1e-50, 1e-50, 0.1802433115377602, 0.09596130634503, 0.04924410173583071, 0.0528738027630882, 0.05563515426229366, 0.05408949210236104, 0.056451402818662566, 0.0616788669775358, 0.05791023002578998, 0.04990404737715025, 0.05156259497572963, 0.044581064770191296, 0.03272809371228107, 0.02647597711030644, 0.02375804308750358, 0.019650749819817473, 0.01354625263761169, 0.011800870086227109, 0.010046804039562004, 0.008188536049530657, 0.005913460286034335, 0.005679005913460286, 0.004480683564748482, 0.003673118503660093, 0.0036991689895016543, 0.0037425864659042557, 0.0026571495558392164, 0.002153506829569038, 0.0019972039145196726, 0.0012243728345533644, 0.0007207301082831862, 0.0010941204053455597, 0.0007120466130026658, 0.0007815145752468283, 0.0007467805941247471, 0.0005470602026727799, 0.000494959230989658, 0.00026918835369612974, 0.0002605048584156094, 0.0002605048584156094, 0.0005383767073922595, 0.00019972039145196724, 6.946796224416251e-05, 0.0001128854386467641, 0.00012156893392728441, 8.683495280520315e-05, 0.00012156893392728441, 2.6050485841560945e-05, 6.0784466963642205e-05, 8.683495280520315e-05, 2.6050485841560945e-05, 5.210097168312189e-05, 5.210097168312189e-05, 6.0784466963642205e-05, 6.0784466963642205e-05, 5.210097168312189e-05, 4.341747640260157e-05, 0.00010420194336624378, 9.551844808572346e-05, 3.473398112208126e-05, 4.341747640260157e-05, 4.341747640260157e-05, 1.736699056104063e-05, 3.473398112208126e-05, 7.815145752468283e-05, 2.6050485841560945e-05, 5.210097168312189e-05, 3.473398112208126e-05, 1.736699056104063e-05, 1.736699056104063e-05, 1.736699056104063e-05, 1e-50, 8.683495280520314e-06, 1.736699056104063e-05, 1e-50, 1e-50, 8.683495280520314e-06, 8.683495280520314e-06, 1.736699056104063e-05, 1e-50, 1e-50, 8.683495280520314e-06, 1.736699056104063e-05, 8.683495280520314e-06, 8.683495280520314e-06, 8.683495280520314e-06, 1e-50, 8.683495280520314e-06, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1.736699056104063e-05]
    prbrun = [ 0.19596998046674205, 0.10458174839793016, 0.10431445118398958, 0.11898838285185566, 0.12412871388917446, 0.11188787224563929, 0.07712552688393133, 0.05459031561632569, 0.036523765463829204, 0.02563997121414619, 0.016791748055241425, 0.01026010075048833, 0.0064356944587231416, 0.004523491312840547, 0.003276104314451184, 0.001528391761762791, 0.001165141701792262, 0.000671669922209657, 0.0004523491312840547, 0.0003632500599705288, 0.0002467358897913026, 0.0001644905931942017, 0.00015763681847777664, 4.7976423014975495e-05, 7.539152188067578e-05, 2.7415098865700285e-05, 1.3707549432850143e-05, 6.853774716425071e-06, 1e-50]
    #prbrun = [0.19596998046674205, 0.10458174839793016, 0.10431445118398958, 0.11898838285185566, 0.12412871388917446, 0.11188787224563929, 0.07712552688393133, 0.05459031561632569, 0.036523765463829204, 0.02563997121414619, 0.016791748055241425, 0.01026010075048833, 0.0064356944587231416, 0.004523491312840547, 0.003276104314451184, 0.001528391761762791, 0.001165141701792262, 0.000671669922209657, 0.0004523491312840547, 0.0003632500599705288, 0.0002467358897913026, 0.0001644905931942017, 0.00015763681847777664, 4.7976423014975495e-05, 7.539152188067578e-05, 2.7415098865700285e-05, 1.3707549432850143e-05, 6.853774716425071e-06, 1e-50, 6.853774716425071e-06, 1e-50, 6.853774716425071e-06, 1e-50, 1.3707549432850143e-05, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 6.853774716425071e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 6.853774716425071e-06]
    prcrun = [ 0.1254604541674868, 0.1787104471376522, 0.1489883063722438, 0.14482664427831718, 0.11945898392778954, 0.07856141464374806, 0.05453946549154612, 0.038881011010729535, 0.02815147485930288, 0.02177641912275698, 0.015019743792656031, 0.011014746584504638, 0.008620585765990865, 0.006266595431009203, 0.00457541807430736, 0.003342184229871575, 0.0026592860098257004, 0.0021210015304954227, 0.0015345124709266126, 0.001004262088302757, 0.001004262088302757, 0.0006587959299266085, 0.0005262333342706447, 0.0003776025452018366, 0.0004499094155596351, 0.0002771763363715609, 0.0003093127231972491, 0.0001606819341284411, 9.640916047706467e-05]
    #prcrun = [0.1254604541674868, 0.1787104471376522, 0.1489883063722438, 0.14482664427831718, 0.11945898392778954, 0.07856141464374806, 0.05453946549154612, 0.038881011010729535, 0.02815147485930288, 0.02177641912275698, 0.015019743792656031, 0.011014746584504638, 0.008620585765990865, 0.006266595431009203, 0.00457541807430736, 0.003342184229871575, 0.0026592860098257004, 0.0021210015304954227, 0.0015345124709266126, 0.001004262088302757, 0.001004262088302757, 0.0006587959299266085, 0.0005262333342706447, 0.0003776025452018366, 0.0004499094155596351, 0.0002771763363715609, 0.0003093127231972491, 0.0001606819341284411, 9.640916047706467e-05, 8.837506377064261e-05, 7.23068703577985e-05, 3.213638682568822e-05, 4.4187531885321305e-05, 6.025572529816541e-05, 4.820458023853233e-05, 5.623867694495439e-05, 3.213638682568822e-05, 8.034096706422055e-06, 2.0085241766055138e-05, 2.8119338472477195e-05, 8.034096706422055e-06, 2.0085241766055138e-05, 8.034096706422055e-06, 8.034096706422055e-06, 4.0170483532110275e-06, 2.4102290119266166e-05]
    for i in range(len(prarun)):
        prarun[i] = log(prarun[i])
    for i in range(len(prbrun)):
        prbrun[i] = log(prbrun[i])
    for i in range(len(prcrun)):
        prcrun[i] = log(prcrun[i])
    neginf = -1000000
    #transition = [.0434,.9565,.0603,.9397,.0922,.1219,.7859,.8727,.6820]
    #transition = [.0343,.0955,.0761,.1186,.9238,.9656,.7859,.8727,.6820]
    #transition = [0.004238388638016391, 0.09914768617452506, 0.01225107054769902, 0.21114151989886346, .07174741019132236, 0.10177111316260486, .8264814766460727, 0.8966139251874585, 0.7766074095534375]
    transition = [.01,.01,.01,.01,.01,.01, 1, 1, 1]
    #old nr30.seqs transition = [0.004418553720818186, 0.09726512198326837, 0.012964664875505504, 0.2021373741460104, 0.0922407556667696, 0.12189608291620897, 0.7858631614170214, 0.8983163242959135, 0.7848979609784841]
    for i in range(len(transition)):
        transition[i] = transcoefficient * log(transition[i])
    matrix = [[[0 for i in range(maxAhistory)],[0 for j in range(maxBhistory)],[0 for k in range(maxChistory)]]for l in range(len(prxati))]
    for k in range(minArun - 1):
        for i in range(len(prxati)):
            matrix[i][0][k] = neginf
    for k in range(minBrun - 1):
        for i in range(len(prxati)):
            matrix[i][1][k] = neginf
    for k in range(minCrun - 1):
        for i in range(len(prxati)):
            matrix[i][2][k] = neginf

    matrix[0][0][:] = [neginf for i in range(maxAhistory)]
    matrix[0][1][:] = [neginf for i in range(maxBhistory)]
    matrix[0][2][:] = [neginf for i in range(maxChistory)]
    matrix[0][2][0] = x * prxati[0][2] + y * prarun[0]
    for i in range(1,len(prxati)):
        for j in range(maxAhistory):
            matrix[i][0][j] = neginf
        for j in range(maxBhistory):
            matrix[i][1][j] = neginf
        for j in range(maxChistory):
            matrix[i][2][j] = neginf

        #next we will do the A ones
        for j in range(minArun,min(i+1,maxAhistory - 1)):
            matrix[i][0][j] = z * matrix[i-1][0][j-1] + x * prxati[i][0] + y * prarun[j] - y * prarun[j-1] + w * transition[7]
        if i >= maxAhistory - 1:
            matrix[i][0][maxAhistory - 1] = x * sum([prxati[j][0] for j in range(i-maxAhistory+minArun,i)]) + max([matrix[i - maxAhistory+minArun][0][l] * z + y * prarun[maxAhistory - 1] - y * prarun[l] for l in range(minArun - 1,maxAhistory)]) + w * transition[7]
        if i >= minArun:
            #this means that we set the value at lambda_A
            matrix[i][0][minArun - 1] = x * sum([prxati[i][0] for i in range(i-minArun + 1,i + 1)]) + z * max(max(matrix[i-minArun][1][minBrun - 1:maxBhistory]) + w * transition[2],max(matrix[i-minArun][2][minCrun - 1:maxChistory]) + w * transition[4]) + y * prarun[minArun - 1]
        if i == minArun - 1:
            matrix[i][0][minArun-1] = x * sum([prxati[j][0] for j in range(minArun)]) + y * prarun[minArun-1]


        #now we go to the B ones
        for j in range(minBrun,min(i+1,maxBhistory - 1)):
            matrix[i][1][j] = z * matrix[i-1][1][j-1] + x * prxati[i][1] + y * prbrun[j] - y * prbrun[j-1] + w * transition[8]
        if i >= maxBhistory - 1:
            matrix[i][1][maxBhistory - 1] = x * sum([prxati[j][1] for j in range(i-maxBhistory+minBrun,i)]) + max([matrix[i - maxBhistory+minBrun][1][l] * z + y * prbrun[maxBhistory - 1] - y * prbrun[l] for l in range(minBrun - 1,maxBhistory)]) + w * transition[8]
        if i >= minBrun:
            matrix[i][1][minBrun - 1] = x * sum([prxati[i][1] for i in range(i-minBrun + 1,i + 1)]) + z * max(max(matrix[i-minBrun][0][minArun - 1:maxAhistory]) + w * transition[0],max(matrix[i-minBrun][2][minCrun - 1:maxChistory]) + w * transition[5]) + y * prbrun[minBrun - 1]
        if i == minBrun - 1:
            matrix[i][1][minBrun-1] = x * sum([prxati[j][1] for j in range(minBrun)]) + y * prbrun[minBrun-1]
        
        #now for C ones
        for j in range(minCrun,min(i+1,maxChistory - 1)):
            matrix[i][2][j] = z * matrix[i-1][2][j-1] + x * prxati[i][2] + y * prcrun[j] - y * prcrun[j-1] + w * transition[6]
        if i >= maxChistory - 1:
            matrix[i][2][maxChistory - 1] = x * sum([prxati[j][2] for j in range(i-maxChistory+minCrun,i)]) + max([matrix[i - maxChistory+minCrun][2][l] * z + y * prcrun[maxChistory - 1] - y * prcrun[l] for l in range(minCrun - 1,maxChistory)]) + w * transition[6]
        if i >= minCrun:
            matrix[i][2][minCrun - 1] = x * sum([prxati[i][2] for i in range(i-minCrun + 1,i + 1)]) + z * max(max(matrix[i-minCrun][0][minArun - 1:maxAhistory]) + w * transition[1],max(matrix[i-minCrun][1][minBrun - 1:maxBhistory]) + w * transition[3]) + y * prcrun[minCrun - 1]
        if i == minCrun - 1:
            matrix[i][2][minCrun-1] = x * sum([prxati[j][2] for j in range(minCrun)]) + y * prcrun[minCrun-1]

    for i in range(maxAhistory):
        matrix[-1][0][i] = neginf
    for i in range(maxBhistory):
        matrix[-1][1][i] = neginf
    return matrix





def fillmatrixlength(transition):
    global prxati
    transcoefficient = 1
    minArun = 3
    minBrun = 1
    minCrun = 1
    import sys
    w = float(sys.argv[4])
    #w = .00001
    x = float(sys.argv[2])
    y = float(sys.argv[3])
    #y=0
    z = 1 
    #z = float(sys.argv[4])
    #prarun = [1e-50, 1e-50, 0.1802433115377602, 0.09596130634503, 0.04924410173583071, 0.0528738027630882, 0.05563515426229366, 0.05408949210236104, 0.056451402818662566, 0.0616788669775358, 0.05791023002578998, 0.04990404737715025, 0.05156259497572963, 0.044581064770191296, 0.03272809371228107, 0.02647597711030644, 0.02375804308750358, 0.019650749819817473, 0.01354625263761169, 0.011800870086227109, 0.010046804039562004, 0.008188536049530657, 0.005913460286034335, 0.005679005913460286, 0.004480683564748482, 0.003673118503660093, 0.0036991689895016543, 0.0037425864659042557, 0.0026571495558392164, 0.002153506829569038, 0.0019972039145196726, 0.0012243728345533644, 0.0007207301082831862, 0.0010941204053455597, 0.0007120466130026658, 0.0007815145752468283, 0.0007467805941247471, 0.0005470602026727799, 0.000494959230989658]
    prarun = [1e-50, 1e-50, 0.1802433115377602, 0.09596130634503, 0.04924410173583071, 0.0528738027630882, 0.05563515426229366, 0.05408949210236104, 0.056451402818662566, 0.0616788669775358, 0.05791023002578998, 0.04990404737715025, 0.05156259497572963, 0.044581064770191296, 0.03272809371228107, 0.02647597711030644, 0.02375804308750358, 0.019650749819817473, 0.01354625263761169, 0.011800870086227109, 0.010046804039562004, 0.008188536049530657, 0.005913460286034335, 0.005679005913460286, 0.004480683564748482, 0.003673118503660093, 0.0036991689895016543, 0.0037425864659042557, 0.0026571495558392164, 0.002153506829569038, 0.0019972039145196726, 0.0012243728345533644, 0.0007207301082831862, 0.0010941204053455597, 0.0007120466130026658, 0.0007815145752468283, 0.0007467805941247471, 0.0005470602026727799, 0.000494959230989658, 0.00026918835369612974, 0.0002605048584156094, 0.0002605048584156094, 0.0005383767073922595, 0.00019972039145196724, 6.946796224416251e-05, 0.0001128854386467641, 0.00012156893392728441, 8.683495280520315e-05, 0.00012156893392728441, 2.6050485841560945e-05, 6.0784466963642205e-05, 8.683495280520315e-05, 2.6050485841560945e-05, 5.210097168312189e-05, 5.210097168312189e-05, 6.0784466963642205e-05, 6.0784466963642205e-05, 5.210097168312189e-05, 4.341747640260157e-05, 0.00010420194336624378, 9.551844808572346e-05, 3.473398112208126e-05, 4.341747640260157e-05, 4.341747640260157e-05, 1.736699056104063e-05, 3.473398112208126e-05, 7.815145752468283e-05, 2.6050485841560945e-05, 5.210097168312189e-05, 3.473398112208126e-05, 1.736699056104063e-05, 1.736699056104063e-05, 1.736699056104063e-05, 1e-50, 8.683495280520314e-06, 1.736699056104063e-05, 1e-50, 1e-50, 8.683495280520314e-06, 8.683495280520314e-06, 1.736699056104063e-05, 1e-50, 1e-50, 8.683495280520314e-06, 1.736699056104063e-05, 8.683495280520314e-06, 8.683495280520314e-06, 8.683495280520314e-06, 1e-50, 8.683495280520314e-06, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1.736699056104063e-05]
    #prbrun = [ 0.19596998046674205, 0.10458174839793016, 0.10431445118398958, 0.11898838285185566, 0.12412871388917446, 0.11188787224563929, 0.07712552688393133, 0.05459031561632569, 0.036523765463829204, 0.02563997121414619, 0.016791748055241425, 0.01026010075048833, 0.0064356944587231416, 0.004523491312840547, 0.003276104314451184, 0.001528391761762791, 0.001165141701792262, 0.000671669922209657, 0.0004523491312840547, 0.0003632500599705288, 0.0002467358897913026, 0.0001644905931942017, 0.00015763681847777664, 4.7976423014975495e-05, 7.539152188067578e-05, 2.7415098865700285e-05, 1.3707549432850143e-05, 6.853774716425071e-06, 1e-50]
    prbrun = [0.19596998046674205, 0.10458174839793016, 0.10431445118398958, 0.11898838285185566, 0.12412871388917446, 0.11188787224563929, 0.07712552688393133, 0.05459031561632569, 0.036523765463829204, 0.02563997121414619, 0.016791748055241425, 0.01026010075048833, 0.0064356944587231416, 0.004523491312840547, 0.003276104314451184, 0.001528391761762791, 0.001165141701792262, 0.000671669922209657, 0.0004523491312840547, 0.0003632500599705288, 0.0002467358897913026, 0.0001644905931942017, 0.00015763681847777664, 4.7976423014975495e-05, 7.539152188067578e-05, 2.7415098865700285e-05, 1.3707549432850143e-05, 6.853774716425071e-06, 1e-50, 6.853774716425071e-06, 1e-50, 6.853774716425071e-06, 1e-50, 1.3707549432850143e-05, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 6.853774716425071e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 6.853774716425071e-06]
    #prcrun = [ 0.1254604541674868, 0.1787104471376522, 0.1489883063722438, 0.14482664427831718, 0.11945898392778954, 0.07856141464374806, 0.05453946549154612, 0.038881011010729535, 0.02815147485930288, 0.02177641912275698, 0.015019743792656031, 0.011014746584504638, 0.008620585765990865, 0.006266595431009203, 0.00457541807430736, 0.003342184229871575, 0.0026592860098257004, 0.0021210015304954227, 0.0015345124709266126, 0.001004262088302757, 0.001004262088302757, 0.0006587959299266085, 0.0005262333342706447, 0.0003776025452018366, 0.0004499094155596351, 0.0002771763363715609, 0.0003093127231972491, 0.0001606819341284411, 9.640916047706467e-05]
    prcrun = [0.1254604541674868, 0.1787104471376522, 0.1489883063722438, 0.14482664427831718, 0.11945898392778954, 0.07856141464374806, 0.05453946549154612, 0.038881011010729535, 0.02815147485930288, 0.02177641912275698, 0.015019743792656031, 0.011014746584504638, 0.008620585765990865, 0.006266595431009203, 0.00457541807430736, 0.003342184229871575, 0.0026592860098257004, 0.0021210015304954227, 0.0015345124709266126, 0.001004262088302757, 0.001004262088302757, 0.0006587959299266085, 0.0005262333342706447, 0.0003776025452018366, 0.0004499094155596351, 0.0002771763363715609, 0.0003093127231972491, 0.0001606819341284411, 9.640916047706467e-05, 8.837506377064261e-05, 7.23068703577985e-05, 3.213638682568822e-05, 4.4187531885321305e-05, 6.025572529816541e-05, 4.820458023853233e-05, 5.623867694495439e-05, 3.213638682568822e-05, 8.034096706422055e-06, 2.0085241766055138e-05, 2.8119338472477195e-05, 8.034096706422055e-06, 2.0085241766055138e-05, 8.034096706422055e-06, 8.034096706422055e-06, 4.0170483532110275e-06, 2.4102290119266166e-05]
    #prarun = [1e-50, 1e-50, 1e-50, 0.18770974269200646, 0.10134270759875404, 0.05219278195737403, 0.05274559902099406, 0.05499348549977379, 0.05348505090523794, 0.05488717452600071, 0.05946327133041096, 0.054520992283004534, 0.0491900875648054, 0.04869160722111383, 0.04270157446552158, 0.03137354959347865, 0.027130560506890722, 0.024261345448059413, 0.01914778760957413, 0.013689309722847291, 0.014118097317065393, 0.009740447663698228, 0.0075823348961046476, 0.00584592232447763, 0.005745517515914162, 0.004593815300039099, 0.0035082621345117315, 0.0029672574013109325, 0.0032661093609174874, 0.002288048402205126, 0.0018580795749451023, 0.001814373952393946, 0.0015545026831708547, 0.0009804234247962078, 0.000841037925849277, 0.0007111022912377314, 0.000532736101907337, 0.000532736101907337, 0.00040516293337963767, 0.0003071205909000169, 0.00023860907446847466, 0.00020435331625270355, 0.00025632757009732177, 0.00019135975279154898, 0.00018663482062385641, 0.0001322981006953919, 0.00012757316852769933, 0.00012521070244385303, 0.0001263919354857762, 0.0001181233041923142, 9.804234247962077e-05, 0.00010867343985692906, 9.213617727000507e-05, 8.623001206038936e-05, 0.00010867343985692906, 0.00011930453723423733, 5.433671992846453e-05, 5.1974253844618246e-05, 5.551795297038767e-05, 5.551795297038767e-05, 5.433671992846453e-05, 5.315548688654139e-05, 3.543699125769426e-05, 3.4255758215771116e-05, 3.3074525173847974e-05, 2.5987126922309123e-05, 3.071205909000169e-05, 3.3074525173847974e-05, 3.543699125769426e-05, 3.071205909000169e-05, 2.2443427796539697e-05, 1.771849562884713e-05, 1.181233041923142e-05, 1.889972867077027e-05, 1.6537262586923987e-05, 1.771849562884713e-05, 9.449864335385136e-06, 1.4174796503077703e-05, 4.724932167692568e-06, 1.181233041923142e-05, 9.449864335385136e-06, 1.181233041923142e-05, 8.268631293461994e-06, 8.268631293461994e-06, 3.543699125769426e-06, 7.087398251538852e-06, 1.0631097377308278e-05, 4.724932167692568e-06, 5.90616520961571e-06, 4.724932167692568e-06, 1.0631097377308278e-05, 8.268631293461994e-06, 1.181233041923142e-05, 3.543699125769426e-06, 5.90616520961571e-06, 5.90616520961571e-06, 1.181233041923142e-06, 1e-50, 3.543699125769426e-06, 9.449864335385136e-06, 3.543699125769426e-06, 2.362466083846284e-06, 3.543699125769426e-06, 3.543699125769426e-06, 2.362466083846284e-06, 1.181233041923142e-06, 3.543699125769426e-06, 7.087398251538852e-06, 2.362466083846284e-06, 8.268631293461994e-06, 7.087398251538852e-06, 1e-50, 1.181233041923142e-06, 1e-50, 1e-50, 1.181233041923142e-06, 2.362466083846284e-06, 1.181233041923142e-06, 1e-50, 2.362466083846284e-06, 4.0161923425386826e-05]
    #prbrun = [1e-50, 1e-50, 1e-50, 0.15947815394074472, 0.17343760007173142, 0.17359899576016088, 0.16720080954027847, 0.10836823835325161, 0.0780040733197556, 0.04829061471262601, 0.03632811999641343, 0.020656086282647403, 0.015796282775493473, 0.006647965261499443, 0.004054105983168735, 0.002798806184272887, 0.0015934621937004445, 0.0010529147292779465, 0.0008018547694987767, 0.0006686392806363601, 0.0003484097401017049, 0.00021263241491501107, 0.0002190370057257042, 0.00018445221534796143, 0.00010119253480895106, 5.251764464768346e-05, 3.970846302629725e-05, 1.6651936107802074e-05, 1.152826345924759e-05, 5.1236726485544835e-06, 1.2809181621386209e-06, 3.842754486415863e-06, 3.842754486415863e-06, 7.685508972831726e-06, 8.966427134970346e-06, 1.2809181621386209e-06, 2.5618363242772417e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1.2809181621386209e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1.2809181621386209e-06]
    #prcrun = [1e-50, 0.09092330629053649, 0.14045688926591943, 0.1221781532608042, 0.12906789398548363, 0.10697552220664158, 0.07632209522061034, 0.057226694184268755, 0.04574970166403414, 0.038967583462849655, 0.032247506481091676, 0.024351178624900688, 0.020785074655212253, 0.015511571003959876, 0.013664895115520118, 0.011552961360348698, 0.009686027114545184, 0.007904557813870011, 0.007065102130609428, 0.005563451391021116, 0.004874034166985841, 0.004924680060395228, 0.003844656383440059, 0.0031634691170838096, 0.00291657038671305, 0.0024588581250257185, 0.002286662087433804, 0.0019574637802727913, 0.0016909397662058945, 0.0014598678775255682, 0.0012465220515385273, 0.0011585248117397182, 0.0010198816785315222, 0.0008989646080166118, 0.0007938743791921347, 0.0007400631174446614, 0.0006318075202820976, 0.0005419110594804365, 0.0004937974607415192, 0.0004608776300254179, 0.0004165624732922046, 0.000390606452919894, 0.0003703480955561394, 0.00033236367549909943, 0.0003197022021467528, 0.0002747539717459222, 0.00023360418335079561, 0.00021461197332227564, 0.00020258357363754634, 0.00015826841690433306, 0.00016903066925382773, 0.00012978010186155312, 0.00015890149057195041, 0.00013864313320819578, 0.00011901784951205847, 0.00011142096550065049, 8.736416613119185e-05, 9.622719747783451e-05, 8.60980187959572e-05, 0.00010129178681877317, 8.166650312263587e-05, 7.913420845216653e-05, 7.723498744931454e-05, 7.217039810837588e-05, 0.0001000256394835385, 5.950892475602923e-05, 6.520658776458522e-05, 5.381126174747325e-05, 4.051671472750927e-05, 4.3682083065595925e-05, 4.051671472750927e-05, 4.1149788395126597e-05, 4.684745140368259e-05, 4.937974607415192e-05, 3.9883641059891935e-05, 4.937974607415192e-05, 3.355290438371861e-05, 3.038753604563195e-05, 3.355290438371861e-05, 3.102060971324928e-05, 2.9121388710397286e-05, 2.7222167707545287e-05, 2.5956020372310623e-05, 2.5956020372310623e-05, 2.3423725701841294e-05, 2.4056799369458628e-05, 2.9121388710397286e-05, 1.7726062693285306e-05, 2.0258357363754635e-05, 1.3294547019963978e-05, 1.6459915358050638e-05, 2.0891431031371965e-05, 1.899221002851997e-05, 1.96252836961373e-05, 1.899221002851997e-05, 1.8359136360902636e-05, 1.7092989025667972e-05, 1.7726062693285306e-05, 1.2028399684729314e-05, 1.899221002851997e-05, 1.392762068758131e-05, 1.2661473352346646e-05, 1.4560694355198643e-05, 1.076225234949465e-05, 1.392762068758131e-05, 1.392762068758131e-05, 9.496105014259985e-06, 1.392762068758131e-05, 8.863031346642653e-06, 7.596884011407988e-06, 5.697663008555991e-06, 1.5193768022815975e-05, 6.963810343790655e-06, 8.229957679025319e-06, 1.076225234949465e-05, 8.229957679025319e-06, 8.863031346642653e-06, 5.697663008555991e-06, 6.963810343790655e-06, 8.863031346642653e-06, 0.0004304900939797881]
    maxAhistory = len(prarun)
    maxBhistory = len(prbrun)
    maxChistory = len(prcrun)
    #print('maxchist',len(prcrun))
    for i in range(len(prarun)):
        prarun[i] = log(prarun[i])
    for i in range(len(prbrun)):
        prbrun[i] = log(prbrun[i])
    for i in range(len(prcrun)):
        prcrun[i] = log(prcrun[i])
    neginf = -1000000000
    #transition = [.0434,.9565,.0603,.9397,.0922,.1219,.7859,.8727,.6820]
    #transition = [.0343,.0955,.0761,.1186,.9238,.9656,.7859,.8727,.6820]
    #transition = [0.004418553720818186, 0.09726512198326837, 0.012964664875505504, 0.2021373741460104, 0.0922407556667696, 0.12189608291620897, 0.7858631614170214, 0.8983163242959135, 0.7848979609784841]
    #transition = gettransitions()
    matrix = [[[0 for i in range(maxAhistory)],[0 for j in range(maxBhistory)],[0 for k in range(maxChistory)]]for l in range(len(prxati))]
    for k in range(minArun - 1):
        for i in range(len(prxati)):
            matrix[i][0][k] = neginf
    for k in range(minBrun - 1):
        for i in range(len(prxati)):
            matrix[i][1][k] = neginf
    for k in range(minCrun - 1):
        for i in range(len(prxati)):
            matrix[i][2][k] = neginf

    #print(len(matrix),': matrix length')
    #print(len(matrix[0]),': matrix[0] length')
    #print(len(matrix[0][0]),': matrix[0][0] length')
    matrix[0][0][:] = [neginf for i in range(maxAhistory)]
    matrix[0][1][:] = [neginf for i in range(maxBhistory)]
    matrix[0][2][:] = [neginf for i in range(maxChistory)]
    matrix[0][2][0] = x * prxati[0][2] + y * prcrun[0]
    matrix[0][1][0] = x * prxati[0][1] + y * prbrun[0]
    for i in range(1,len(prxati)):
        for j in range(maxAhistory):
            matrix[i][0][j] = neginf
        for j in range(maxBhistory):
            matrix[i][1][j] = neginf
        for j in range(maxChistory):
            matrix[i][2][j] = neginf

        #next we will do the A ones
        for j in range(minArun,min(i+1,maxAhistory - 1)):
            matrix[i][0][j] = z * matrix[i-1][0][j-1] + x * prxati[i][0] + y * prarun[j] - y * prarun[j-1] + w * transition[0][j]
        if i >= maxAhistory - 1:
            matrix[i][0][maxAhistory - 1] = x * sum([prxati[j][0] for j in range(i-maxAhistory+minArun,i)]) + max([matrix[i - maxAhistory+minArun][0][l] * z + y * prarun[maxAhistory - 1] - y * prarun[l] for l in range(minArun - 1,maxAhistory)]) + w * sum([transition[0][m] for m in range(minArun,maxAhistory)])
        if i >= minArun:
            #this means that we set the value at lambda_A
            matrix[i][0][minArun - 1] = x * sum([prxati[p][0] for p in range(i-minArun + 1,i + 1)]) + z * max(max([matrix[i-minArun][1][l] + w * transition[3][l] for l in range(minBrun - 1,maxBhistory)]),max([matrix[i-minArun][2][m] + w * transition[6][m] for m in range(minCrun - 1,maxChistory)])) + y * prarun[minArun - 1] + w * sum([transition[0][m] for m in range(minArun - 1)])
        if i == minArun - 1:
            matrix[i][0][minArun-1] = x * sum([prxati[j][0] for j in range(minArun)]) + y * prarun[minArun-1]


        #now we go to the B ones
        for j in range(minBrun,min(i+1,maxBhistory - 1)):
            matrix[i][1][j] = z * matrix[i-1][1][j-1] + x * prxati[i][1] + y * prbrun[j] - y * prbrun[j-1] + w * transition[4][j]
        if i >= maxBhistory - 1:
            matrix[i][1][maxBhistory - 1] = x * sum([prxati[j][1] for j in range(i-maxBhistory+minBrun,i)]) + max([matrix[i - maxBhistory+minBrun][1][l] * z + y * prbrun[maxBhistory - 1] - y * prbrun[l] for l in range(minBrun - 1,maxBhistory)]) + w * sum([transition[4][m] for m in range(minBrun,maxBhistory)])
        if i >= minBrun:
            matrix[i][1][minBrun - 1] = x * sum([prxati[l][1] for l in range(i-minBrun + 1,i + 1)]) + z * max(max([matrix[i-minBrun][0][l] + w * transition[1][l] for l in range(minArun-1,maxAhistory)]),max([matrix[i-minBrun][2][m] + w * transition[7][m] for m in range(minCrun-1,maxChistory)])) + y * prbrun[minBrun - 1] + w * sum([transition[4][m] for m in range(minBrun - 1)])
        if i == minBrun - 1:
            matrix[i][1][minBrun-1] = x * sum([prxati[j][1] for j in range(minBrun)]) + y * prbrun[minBrun-1]
        
        #now for C ones
        for j in range(minCrun,min(i+1,maxChistory)):
            matrix[i][2][j] = z * matrix[i-1][2][j-1] + x * prxati[i][2] + y * prcrun[j] - y * prcrun[j-1] + w * transition[8][j]
        if i >= maxChistory:
            matrix[i][2][maxChistory - 1] = x * sum([prxati[j][2] for j in range(i-maxChistory+minCrun,i-1)]) + max([matrix[i - maxChistory+minCrun][2][l] * z + y * prcrun[maxChistory - 1] - y * prcrun[l] for l in range(minCrun - 1,maxChistory)]) + w * sum([transition[8][m] for m in range(minCrun,maxChistory)])
        if i >= minCrun:
            matrix[i][2][minCrun - 1] = x * sum([prxati[l][2] for l in range(i-minCrun + 1,i + 1)]) + z * max(max([matrix[i-minCrun][0][l] + w * transition[2][l] for l in range(minArun-1,maxAhistory)]),max([matrix[i-minCrun][1][m] + w * transition[5][m] for m in range(minBrun-1,maxBhistory)])) + y * prcrun[minCrun - 1] + w * sum([transition[8][m] for m in range(minCrun-1 )])
        if i == minCrun - 1:
            matrix[i][2][minCrun-1] = x * sum([prxati[j][2] for j in range(minCrun)]) + y * prcrun[minCrun-1]

    #for i in range(maxAhistory):
        #matrix[-1][0][i] = neginf
    #for i in range(maxBhistory):
        #matrix[-1][1][i] = neginf
    return matrix





def getmax(k,minrunlens,maxhistories,transition,matrix,w,i):
    maxes = []
    for j in range(len(minrunlens)):
        if j != k:
            maxes.append(max(matrix[i-minrunlens[k]][j][l] + w * transition[k][j][l] for l in range(minrunlens[j] - 1,maxhistories[j])))
    return max(maxes)




def fillmatrixlength8(transition):
    global prxati
    transcoefficient = 1
    #             G,H,I,E,B,L,T,S
    minrunlens = [3,4,5,2,1,1,1,1]
    import sys
    w = float(sys.argv[4])
    #w = .00001
    x = float(sys.argv[2])
    y = float(sys.argv[3])
    #y=0
    z = 1 
    #z = float(sys.argv[4])
    prruns = [[1e-50, 2.036991770553247e-05, 0.7961500855536544, 0.12315652244764931, 0.03791249083353703, 0.03010673836877699, 0.008384258127597164, 0.0018332925934979223, 0.0012099731117086287, 0.0009084983296667481, 0.0003177707162063065],
    [1e-50, 1e-50, 1e-50, 0.10458138738065767, 0.059679654775391014, 0.06441785018659538, 0.06983050533548585, 0.0710574973067912, 0.07119948883051459, 0.07669857360254104, 0.0711392967715449, 0.0631569037204866, 0.06140361707972824, 0.05374842188768471, 0.0401203223824943, 0.03362575355827672, 0.029285751767948807, 0.023701780758913826, 0.017333769597145353, 0.01747730450699617, 0.012323938227513636, 0.009372983951871047, 0.007085685711022555, 0.006647364050832966, 0.005543842969721851, 0.004349262107092477, 0.003748884903522933, 0.00418103301920281, 0.0028321135438306227, 0.002253343746045073, 0.002170000895163954, 0.0016730305621320953, 0.0012177316545407962, 0.0010248083886122798, 0.0008812734787614635, 0.0006327883122455342, 0.0006497655596472437, 0.0004877100162672898, 0.0003935634624941737, 0.0003055904532307702, 0.0002531153248982137, 0.00029787352259362954, 0.00023305130524164797, 0.00023613807749650422, 0.00017440263239937895, 0.00016205554337995389, 0.00015742538499766948, 0.00015588199887024134, 0.0001342745930862475, 0.00012655766244910684, 0.00017440263239937895, 0.0001960100381833728, 0.000104950256665113, 7.871269249883474e-05, 7.253914798912221e-05, 6.945237573426595e-05, 5.864867284226903e-05, 6.945237573426595e-05, 7.099576186169409e-05, 0.0009414655377311602],
    [1e-50, 1e-50, 1e-50, 1e-50, 0.7898477157360406, 0.11573604060913706, 0.05685279187817259, 0.026395939086294416, 0.007106598984771574, 0.0020304568527918783, 0.0010152284263959391, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 0.0010152284263959391, 2.9790000000001813e-47],
    [3.857803565051386e-05, 0.1441286434199098, 0.13645271655504754, 0.14947666139066104, 0.14843946334645722, 0.14362161780850305, 0.09239659984216073, 0.06645011529321511, 0.041463672717172294, 0.031233879892246037, 0.020177414874808765, 0.010212157151486026, 0.005512250179663423, 0.0034610009126461007, 0.0023786114552516832, 0.0013623557732581465, 0.0009302817739723913, 0.000654724376468721, 0.0005720571572176198, 0.00029429530053392, 0.0007429027436698952],
    [0.9795000980495678, 0.02021329552139743, 0.0002866064290347397],
    [0.5030707938640901, 0.26076971583865494, 0.10817891901517804, 0.04849824028327525, 0.022673289728261197, 0.01259588972155096, 0.008004538873351152, 0.005610360606433397, 0.004273474861988484, 0.003390390423975162, 0.002829053262938929, 0.0023743916696705013, 0.001895642044970634, 0.001592390934985543, 0.0013510804772527256, 0.0011566556521275325, 0.0010078260293547075, 0.0008688897052055096, 0.0007665155716218902, 0.0007871624557059815, 0.000671884019569805, 0.0006241381001253438, 0.0006151050883385538, 0.0005415505637889785, 0.00048262091546563457, 0.0004374558565316848, 0.000367342479329458, 0.0003535778899400638, 0.0003023908231482541, 0.00026926977993002425, 0.0002408803143143987, 0.00021636213946454026, 0.00019614539879886752, 0.00018152052257263616, 0.00019399468170677468, 0.0001531310569570106, 0.0001449583320070578, 0.00013463488996501213, 0.00011957987032036222, 0.00011484829271775796, 0.00010108370332836374, 8.989997444948095e-05, 9.033011786789952e-05, 8.258753633636527e-05, 6.0650221997018245e-05, 6.581194301804108e-05, 6.839280352855248e-05, 6.237079567069253e-05, 4.817606286287974e-05, 5.0756923373391156e-05, 3.3551186636648395e-05, 4.4734915515531186e-05, 3.183061296297411e-05, 3.1400469544555544e-05, 3.914305107608979e-05, 3.699233398399695e-05, 3.527176031032267e-05, 0.0009704035519522944],
    [0.17777273243273925, 0.6452034421957269, 0.11021592534123456, 0.047025156580377654, 0.011071084102181303, 0.006099059905214663, 0.001555799292907597, 0.0006198696395364126, 0.000436930510081661],
    [0.5969703015815335, 0.279153905736141, 0.08772136120975518, 0.026240544520316718, 0.007151684442829388, 0.002006598637030204, 0.0007556038723939822],
    ]
    prarun = [1e-50, 1e-50, 0.1802433115377602, 0.09596130634503, 0.04924410173583071, 0.0528738027630882, 0.05563515426229366, 0.05408949210236104, 0.056451402818662566, 0.0616788669775358, 0.05791023002578998, 0.04990404737715025, 0.05156259497572963, 0.044581064770191296, 0.03272809371228107, 0.02647597711030644, 0.02375804308750358, 0.019650749819817473, 0.01354625263761169, 0.011800870086227109, 0.010046804039562004, 0.008188536049530657, 0.005913460286034335, 0.005679005913460286, 0.004480683564748482, 0.003673118503660093, 0.0036991689895016543, 0.0037425864659042557, 0.0026571495558392164, 0.002153506829569038, 0.0019972039145196726, 0.0012243728345533644, 0.0007207301082831862, 0.0010941204053455597, 0.0007120466130026658, 0.0007815145752468283, 0.0007467805941247471, 0.0005470602026727799, 0.000494959230989658, 0.00026918835369612974, 0.0002605048584156094, 0.0002605048584156094, 0.0005383767073922595, 0.00019972039145196724, 6.946796224416251e-05, 0.0001128854386467641, 0.00012156893392728441, 8.683495280520315e-05, 0.00012156893392728441, 2.6050485841560945e-05, 6.0784466963642205e-05, 8.683495280520315e-05, 2.6050485841560945e-05, 5.210097168312189e-05, 5.210097168312189e-05, 6.0784466963642205e-05, 6.0784466963642205e-05, 5.210097168312189e-05, 4.341747640260157e-05, 0.00010420194336624378, 9.551844808572346e-05, 3.473398112208126e-05, 4.341747640260157e-05, 4.341747640260157e-05, 1.736699056104063e-05, 3.473398112208126e-05, 7.815145752468283e-05, 2.6050485841560945e-05, 5.210097168312189e-05, 3.473398112208126e-05, 1.736699056104063e-05, 1.736699056104063e-05, 1.736699056104063e-05, 1e-50, 8.683495280520314e-06, 1.736699056104063e-05, 1e-50, 1e-50, 8.683495280520314e-06, 8.683495280520314e-06, 1.736699056104063e-05, 1e-50, 1e-50, 8.683495280520314e-06, 1.736699056104063e-05, 8.683495280520314e-06, 8.683495280520314e-06, 8.683495280520314e-06, 1e-50, 8.683495280520314e-06, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1.736699056104063e-05]
    prbrun = [0.19596998046674205, 0.10458174839793016, 0.10431445118398958, 0.11898838285185566, 0.12412871388917446, 0.11188787224563929, 0.07712552688393133, 0.05459031561632569, 0.036523765463829204, 0.02563997121414619, 0.016791748055241425, 0.01026010075048833, 0.0064356944587231416, 0.004523491312840547, 0.003276104314451184, 0.001528391761762791, 0.001165141701792262, 0.000671669922209657, 0.0004523491312840547, 0.0003632500599705288, 0.0002467358897913026, 0.0001644905931942017, 0.00015763681847777664, 4.7976423014975495e-05, 7.539152188067578e-05, 2.7415098865700285e-05, 1.3707549432850143e-05, 6.853774716425071e-06, 1e-50, 6.853774716425071e-06, 1e-50, 6.853774716425071e-06, 1e-50, 1.3707549432850143e-05, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 6.853774716425071e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 6.853774716425071e-06]
    prcrun = [0.1254604541674868, 0.1787104471376522, 0.1489883063722438, 0.14482664427831718, 0.11945898392778954, 0.07856141464374806, 0.05453946549154612, 0.038881011010729535, 0.02815147485930288, 0.02177641912275698, 0.015019743792656031, 0.011014746584504638, 0.008620585765990865, 0.006266595431009203, 0.00457541807430736, 0.003342184229871575, 0.0026592860098257004, 0.0021210015304954227, 0.0015345124709266126, 0.001004262088302757, 0.001004262088302757, 0.0006587959299266085, 0.0005262333342706447, 0.0003776025452018366, 0.0004499094155596351, 0.0002771763363715609, 0.0003093127231972491, 0.0001606819341284411, 9.640916047706467e-05, 8.837506377064261e-05, 7.23068703577985e-05, 3.213638682568822e-05, 4.4187531885321305e-05, 6.025572529816541e-05, 4.820458023853233e-05, 5.623867694495439e-05, 3.213638682568822e-05, 8.034096706422055e-06, 2.0085241766055138e-05, 2.8119338472477195e-05, 8.034096706422055e-06, 2.0085241766055138e-05, 8.034096706422055e-06, 8.034096706422055e-06, 4.0170483532110275e-06, 2.4102290119266166e-05]
    for i in range(len(prarun)):
        prarun[i] = log(prarun[i])
    for i in range(len(prbrun)):
        prbrun[i] = log(prbrun[i])
    for i in range(len(prcrun)):
        prcrun[i] = log(prcrun[i])
    #prruns = [prarun,prarun,prarun,prbrun,prbrun,prcrun,prcrun,prcrun]
    prruns = [[log(thing) for thing in prrun] for prrun in prruns]
    maxhistories = [len(prruns[x]) for x in range(len(prruns))]
    neginf = -1000000000
    #transition = gettransitions()
    matrix = [[[0 for i in range(maxhistories[j])] for j in range(len(maxhistories))] for l in range(len(prxati))]
    for i,minlen in enumerate(minrunlens):
        for k in range(minlen - 1):
            for j in range(len(prxati)):
                matrix[j][i][k] = neginf

    for k,val in enumerate(maxhistories):
        matrix[0][k][:] = [neginf for i in range(val)]
    for i,minlen in enumerate(minrunlens):
        if minlen == 1:
            matrix[0][i][0] = x * prxati[0][i] + y * prruns[i][0]

    for i in range(1,len(prxati)):
        for j,maxhis in enumerate(maxhistories):
            for k in range(maxhis):
                matrix[i][j][k] = neginf

##############################################################################################################################3
        for k,minlen in enumerate(minrunlens):
            maxhis = maxhistories[k]
            prrun = prruns[k]

            #general case where we are continuing a run
            for j in range(minlen,min(i+1,maxhis - 1)):
                matrix[i][k][j] = z * matrix[i-1][k][j-1] + x * prxati[i][k] + y * prrun[j] - y * prrun[j-1] + w * transition[k][k][j]

            #special case where we have met our maximum history
            if i >= maxhis - 1:
                #note that in C, it had that the first range should end at i-1 instead of i
                matrix[i][k][maxhis - 1] = x * sum([prxati[a][k] for a in range(i-maxhis+minlen,i)]) + max([matrix[i - maxhis+minlen][k][l] * z + y * prrun[maxhis - 1] - y * prrun[l] for l in range(minlen - 1,maxhis)]) + w * sum([transition[k][k][m] for m in range(minlen,maxhis)])

            #special case where we are at the min len so we need to look at other ss types and transition
            if i >= minlen:
                matrix[i][k][minlen - 1] = x * sum([prxati[p][k] for p in range(i-minlen + 1,i + 1)]) + z * getmax(k,minrunlens,maxhistories,transition,matrix,w,i) + y * prrun[minlen - 1] + w * sum([transition[k][k][m] for m in range(minlen - 1)])

            #special case where we are filling the matrix at the length of our minlen (here we don't need to transition)
            if i == minlen - 1:
                matrix[i][k][minlen-1] = x * sum([prxati[j][k] for j in range(minlen)]) + y * prrun[minlen-1]

##############################################################################################################################3


    return matrix



def backtrack2(matrix):
    maxAhistory = len(matrix[0][0])
    maxBhistory = len(matrix[0][1])
    maxChistory = len(matrix[0][2])
    minArun = 3
    minBrun = 1
    minCrun = 1
    letteranswer = [' ' for i in range(len(matrix))]
    numberanswer = [0 for i in range(len(matrix))]
    lastA = max(matrix[-1][0][:])
    lastB = max(matrix[-1][1][:])
    lastC = max(matrix[-1][2][:])

    if lastA >= lastB and lastA >= lastC:
        letter = 'A'
        number = matrix[-1][0][:].index(lastA)
    elif lastB >= lastC:
        letter = 'B'
        number = matrix[-1][1][:].index(lastB)
    else:
        letter = 'C'
        number = matrix[-1][2][:].index(lastC)
    letteranswer[-1] = letter
    numberanswer[-1] = number

    rangeiterator = number
    rangeletter = letteranswer[-1]
    maxhistoryhit = False

    for i in range(len(matrix) - 2,-1,-1):
        if numberanswer[i+1] > 0:
            if letteranswer[i+1] == 'A' and numberanswer[i+1] == maxAhistory - 1:
                maxhistoryhit = True
            elif letteranswer[i+1] == 'B' and numberanswer[i+1] == maxBhistory - 1:
                maxhistoryhit = True
            elif letteranswer[i+1] == 'C' and numberanswer[i+1] == maxChistory - 1:
                maxhistoryhit = True
            if maxhistoryhit == True:
                if numberanswer[i+1] == minArun - 1 and letteranswer[i+1] == 'A':
                    letteranswer[i] = 'A'
                    maxval = max(matrix[i][0][:])
                    numberanswer[i] = matrix[i][0][:].index(maxval)
                    maxhistoryhit = False
                elif numberanswer[i+1] == minBrun - 1 and letteranswer[i+1] == 'B':
                    letteranswer[i] = 'B'
                    maxval = max(matrix[i][1][:])
                    numberanswer[i] = matrix[i][1][:].index(maxval)
                    maxhistoryhit = False
                elif numberanswer[i+1] == minCrun - 1 and letteranswer[i+1] == 'C':
                    letteranswer[i] = 'C'
                    maxval = max(matrix[i][2][:])
                    numberanswer[i] = matrix[i][2][:].index(maxval)
                    maxhistoryhit = False
                else:
                    letteranswer[i] = letteranswer[i+1]
                    numberanswer[i] = numberanswer[i+1] - 1
            else:
                letteranswer[i] = letteranswer[i+1]
                numberanswer[i] = numberanswer[i+1] - 1
        else:
            if letteranswer[i+1] == 'A':
                betavalue = max(matrix[i][1][:])
                coilvalue = max(matrix[i][2][:])
                if betavalue >= coilvalue:
                    letteranswer[i] = 'B'
                    numberanswer[i] = matrix[i][1][:].index(betavalue)
                else:
                    letteranswer[i] = 'C'
                    numberanswer[i] = matrix[i][2][:].index(coilvalue)
            elif letteranswer[i+1] == 'B':
                alphavalue = max(matrix[i][0][:])
                coilvalue = max(matrix[i][2][:])
                if alphavalue >= coilvalue:
                    letteranswer[i] = 'A'
                    numberanswer[i] = matrix[i][0][:].index(alphavalue)
                else:
                    letteranswer[i] = 'C'
                    numberanswer[i] = matrix[i][2][:].index(coilvalue)
            elif letteranswer[i+1] == 'C':
                alphavalue = max(matrix[i][0][:])
                betavalue = max(matrix[i][1][:])
                if alphavalue >= betavalue:
                    letteranswer[i] = 'A'
                    numberanswer[i] = matrix[i][0][:].index(alphavalue)
                else:
                    letteranswer[i] = 'B'
                    numberanswer[i] = matrix[i][1][:].index(betavalue)
                

    #print(numberanswer)
    return ''.join(letteranswer)

def backtrack8(matrix):
    numclasses = len(matrix[0])
    maxhistories = [len(matrix[0][k]) for k in range(numclasses)]
    minrunlens = [3,5,7,1,2,1,1,1]
    letteranswer = [' ' for i in range(len(matrix))]
    numberanswer = [0 for i in range(len(matrix))]
    lasts = [max(matrix[-1][k][:]) for k in range(numclasses)]
    letters = ['G','H','I','E','B','L','T','S']
    neginf = -1000000

    maxlast = max(lasts)
    maxindex = lasts.index(maxlast)
    letter = letters[maxindex]
    number = matrix[-1][maxindex][:].index(maxlast)

    letteranswer[-1] = letter
    numberanswer[-1] = number

    rangeiterator = number
    rangeletter = letteranswer[-1]
    maxhistoryhit = False

    for i in range(len(matrix) - 2,-1,-1):
        letternum = letters.index(letteranswer[i+1])
        if numberanswer[i+1] > 0:
            if numberanswer[i+1] == maxhistories[letternum] - 1:
                maxhistoryhit = True

            if maxhistoryhit == True:
                if numberanswer[i+1] == minrunlens[letternum] - 1:
                    letteranswer[i] = letteranswer[i+1]
                    maxval = max(matrix[i][letternum][:])
                    numberanswer[i] = matrix[i][letternum][:].index(maxval)
                    maxhistoryhit = False
                else:
                    letteranswer[i] = letteranswer[i+1]
                    numberanswer[i] = numberanswer[i+1] - 1
            else:
                letteranswer[i] = letteranswer[i+1]
                numberanswer[i] = numberanswer[i+1] - 1
        else:
            othervals = []
            for j in range(8):
                if j == letternum:
                    othervals.append(neginf)
                else:
                    othervals.append(max(matrix[i][j][:]))
            maxothervalue = max(othervals)
            letteranswer[i] = letters[othervals.index(maxothervalue)]
            numberanswer[i] = matrix[i][othervals.index(maxothervalue)][:].index(maxothervalue)

    #print(numberanswer)
    return ''.join(letteranswer)

def backtrack(matrix):
    maxAhistory = 39
    maxBhistory = 28
    maxChistory = 29
    minArun = 3
    minBrun = 3
    minCrun = 1
    letteranswer = [' ' for i in range(len(matrix))]
    numberanswer = [' ' for i in range(len(matrix))]
    lastA = max(matrix[-1][0][:])
    lastB = max(matrix[-1][1][:])
    lastC = max(matrix[-1][2][:])

    if lastA >= lastB and lastA >= lastC:
        letter = 'A'
        number = matrix[-1][0][:].index(lastA)
    elif lastB >= lastC:
        letter = 'B'
        number = matrix[-1][1][:].index(lastB)
    else:
        letter = 'C'
        number = matrix[-1][2][:].index(lastC)
    letteranswer[-1] = letter
    numberanswer[-1] = number

    rangeiterator = 0
    rangeletter = letteranswer[-1]

    for i in range(len(matrix) - 2,-1,-1):
        if rangeiterator > 0:
            letteranswer[i] = letteranswer[i+1]
            numberanswer[i] = numberanswer[i+1]
            rangeiterator -= 1
            #print('im in rangei giving numberanswer {0} '.format(numberanswer[i]))
            continue
        if rangeletter == 'A':
            if numberanswer[i+1] in range(minArun,maxAhistory):
                number = numberanswer[i+1] - 1
                letter = 'A'
                rangeletter = 'A'
            elif numberanswer[i+1] == minArun - 1:
                rangeiterator = minArun - 1
                letter = 'A'
                betavalue = max(matrix[i-minArun][1][:])
                coilvalue = max(matrix[i-minArun][2][:])
                if betavalue >= coilvalue:
                    rangeletter = 'B'
                    number = matrix[i-minArun][1][:].index(betavalue) + 1
                else:
                    rangeletter = 'C'
                    number = matrix[i-minArun][2][:].index(coilvalue) + 1
            #need to do something here to skip back and assign everything to 'A' and move on.
        elif rangeletter == 'B':
            if numberanswer[i+1] in range(minBrun,maxBhistory):
                number = numberanswer[i+1] - 1
                letter = 'B'
                rangeletter = 'B'
            elif numberanswer[i+1] == minBrun - 1:
                rangeiterator = minBrun - 1
                letter = 'B'
                alphavalue = max(matrix[i-minBrun][0][:])
                coilvalue = max(matrix[i-minBrun][2][:])
                if alphavalue >= coilvalue:
                    rangeletter = 'A'
                    number = matrix[i-minBrun][0][:].index(alphavalue) + 1
                else:
                    rangeletter = 'C'
                    number = matrix[i-minBrun][2][:].index(coilvalue) + 1
        elif rangeletter == 'C':
            if numberanswer[i+1] in range(minCrun,maxChistory):
                rangeiterator = minCrun - 1
                number = numberanswer[i+1] - 1
                letter = 'C'
                rangeletter = 'C'
            elif numberanswer[i+1] == minCrun - 1:
                letter = 'C'
                alphavalue = max(matrix[i-minCrun][0][:])
                betavalue = max(matrix[i-minCrun][1][:])
                if alphavalue >= betavalue:
                    rangeletter = 'A'
                    number = matrix[i-minCrun][0][:].index(alphavalue) + 1
                else:
                    rangeletter = 'B'
                    number = matrix[i-minCrun][1][:].index(betavalue) + 1
            
        letteranswer[i] = letter
        numberanswer[i] = number

    #print(numberanswer)
    return ''.join(letteranswer)

import sys
if sys.argv[1] == 'h':
    print('usage: python newdponprotein.py infile outfile type')
    exit()
out = open(sys.argv[2],'w')
import traceback
rand = sys.argv[1]
overallanswer = []

if sys.argv[3] == '8':
    transitions = gettransitions8()
else:
    transitions = gettransitions()
try:
    prxati = getprxati(i,rand)
    if sys.argv[3] == '8':
        matrix = fillmatrixlength8(transitions)
        answer = backtrack8(matrix)
    else:
        matrix = fillmatrixlength(transitions)
        answer = backtrack(matrix)
except:
    print(traceback.format_exc())
    a,b,c = sys.exc_info()
    #print(c.tb_lineno)
    answer = ''

out.write(answer + '\n')
overallanswer.append(answer)


