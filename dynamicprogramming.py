
from numpy import inf
from math import log
import pickle

def getprxati(rand):
    prxati = []
    neginf = -1000000
    fin = open(rand,'r').readlines()
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
        

def gettransitions():
    with open('transitions.txt','rb') as transitionfile:
        transitions = pickle.load(transitionfile)
    return transitions
 
def gettransitions8():
    with open('transitions8.txt','rb') as transitionfile:
        transitions = pickle.load(transitionfile)
    return transitions




def fillmatrixlength(transition):
    global prxati
    transcoefficient = 1
    minArun = 3
    minBrun = 1
    minCrun = 1
    import sys
    w = .1
    x = 1
    y = .1
    z = 1 
    prarun = [1e-50, 1e-50, 0.1802433115377602, 0.09596130634503, 0.04924410173583071, 0.0528738027630882, 0.05563515426229366, 0.05408949210236104, 0.056451402818662566, 0.0616788669775358, 0.05791023002578998, 0.04990404737715025, 0.05156259497572963, 0.044581064770191296, 0.03272809371228107, 0.02647597711030644, 0.02375804308750358, 0.019650749819817473, 0.01354625263761169, 0.011800870086227109, 0.010046804039562004, 0.008188536049530657, 0.005913460286034335, 0.005679005913460286, 0.004480683564748482, 0.003673118503660093, 0.0036991689895016543, 0.0037425864659042557, 0.0026571495558392164, 0.002153506829569038, 0.0019972039145196726, 0.0012243728345533644, 0.0007207301082831862, 0.0010941204053455597, 0.0007120466130026658, 0.0007815145752468283, 0.0007467805941247471, 0.0005470602026727799, 0.000494959230989658, 0.00026918835369612974, 0.0002605048584156094, 0.0002605048584156094, 0.0005383767073922595, 0.00019972039145196724, 6.946796224416251e-05, 0.0001128854386467641, 0.00012156893392728441, 8.683495280520315e-05, 0.00012156893392728441, 2.6050485841560945e-05, 6.0784466963642205e-05, 8.683495280520315e-05, 2.6050485841560945e-05, 5.210097168312189e-05, 5.210097168312189e-05, 6.0784466963642205e-05, 6.0784466963642205e-05, 5.210097168312189e-05, 4.341747640260157e-05, 0.00010420194336624378, 9.551844808572346e-05, 3.473398112208126e-05, 4.341747640260157e-05, 4.341747640260157e-05, 1.736699056104063e-05, 3.473398112208126e-05, 7.815145752468283e-05, 2.6050485841560945e-05, 5.210097168312189e-05, 3.473398112208126e-05, 1.736699056104063e-05, 1.736699056104063e-05, 1.736699056104063e-05, 1e-50, 8.683495280520314e-06, 1.736699056104063e-05, 1e-50, 1e-50, 8.683495280520314e-06, 8.683495280520314e-06, 1.736699056104063e-05, 1e-50, 1e-50, 8.683495280520314e-06, 1.736699056104063e-05, 8.683495280520314e-06, 8.683495280520314e-06, 8.683495280520314e-06, 1e-50, 8.683495280520314e-06, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 8.683495280520314e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1.736699056104063e-05]
    prbrun = [0.19596998046674205, 0.10458174839793016, 0.10431445118398958, 0.11898838285185566, 0.12412871388917446, 0.11188787224563929, 0.07712552688393133, 0.05459031561632569, 0.036523765463829204, 0.02563997121414619, 0.016791748055241425, 0.01026010075048833, 0.0064356944587231416, 0.004523491312840547, 0.003276104314451184, 0.001528391761762791, 0.001165141701792262, 0.000671669922209657, 0.0004523491312840547, 0.0003632500599705288, 0.0002467358897913026, 0.0001644905931942017, 0.00015763681847777664, 4.7976423014975495e-05, 7.539152188067578e-05, 2.7415098865700285e-05, 1.3707549432850143e-05, 6.853774716425071e-06, 1e-50, 6.853774716425071e-06, 1e-50, 6.853774716425071e-06, 1e-50, 1.3707549432850143e-05, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 6.853774716425071e-06, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 1e-50, 6.853774716425071e-06]
    prcrun = [0.1254604541674868, 0.1787104471376522, 0.1489883063722438, 0.14482664427831718, 0.11945898392778954, 0.07856141464374806, 0.05453946549154612, 0.038881011010729535, 0.02815147485930288, 0.02177641912275698, 0.015019743792656031, 0.011014746584504638, 0.008620585765990865, 0.006266595431009203, 0.00457541807430736, 0.003342184229871575, 0.0026592860098257004, 0.0021210015304954227, 0.0015345124709266126, 0.001004262088302757, 0.001004262088302757, 0.0006587959299266085, 0.0005262333342706447, 0.0003776025452018366, 0.0004499094155596351, 0.0002771763363715609, 0.0003093127231972491, 0.0001606819341284411, 9.640916047706467e-05, 8.837506377064261e-05, 7.23068703577985e-05, 3.213638682568822e-05, 4.4187531885321305e-05, 6.025572529816541e-05, 4.820458023853233e-05, 5.623867694495439e-05, 3.213638682568822e-05, 8.034096706422055e-06, 2.0085241766055138e-05, 2.8119338472477195e-05, 8.034096706422055e-06, 2.0085241766055138e-05, 8.034096706422055e-06, 8.034096706422055e-06, 4.0170483532110275e-06, 2.4102290119266166e-05]
    maxAhistory = len(prarun)
    maxBhistory = len(prbrun)
    maxChistory = len(prcrun)
    for i in range(len(prarun)):
        prarun[i] = log(prarun[i])
    for i in range(len(prbrun)):
        prbrun[i] = log(prbrun[i])
    for i in range(len(prcrun)):
        prcrun[i] = log(prcrun[i])
    neginf = -1000000000

    #setting initial boundary conditions
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
    matrix[0][2][0] = x * prxati[0][2] + y * prcrun[0]
    matrix[0][1][0] = x * prxati[0][1] + y * prbrun[0]
    for i in range(1,len(prxati)):
        for j in range(maxAhistory):
            matrix[i][0][j] = neginf
        for j in range(maxBhistory):
            matrix[i][1][j] = neginf
        for j in range(maxChistory):
            matrix[i][2][j] = neginf

        #next we will do the Alpha recurrences
        for j in range(minArun,min(i+1,maxAhistory - 1)):
            matrix[i][0][j] = z * matrix[i-1][0][j-1] + x * prxati[i][0] + y * prarun[j] - y * prarun[j-1] + w * transition[0][j]
        if i >= maxAhistory - 1:
            matrix[i][0][maxAhistory - 1] = x * sum([prxati[j][0] for j in range(i-maxAhistory+minArun,i)]) + max([matrix[i - maxAhistory+minArun][0][l] * z + y * prarun[maxAhistory - 1] - y * prarun[l] for l in range(minArun - 1,maxAhistory)]) + w * sum([transition[0][m] for m in range(minArun,maxAhistory)])
        if i >= minArun:
            #this means that we set the value at lambda_A
            matrix[i][0][minArun - 1] = x * sum([prxati[p][0] for p in range(i-minArun + 1,i + 1)]) + z * max(max([matrix[i-minArun][1][l] + w * transition[3][l] for l in range(minBrun - 1,maxBhistory)]),max([matrix[i-minArun][2][m] + w * transition[6][m] for m in range(minCrun - 1,maxChistory)])) + y * prarun[minArun - 1] + w * sum([transition[0][m] for m in range(minArun - 1)])
        if i == minArun - 1:
            matrix[i][0][minArun-1] = x * sum([prxati[j][0] for j in range(minArun)]) + y * prarun[minArun-1]


        #now we go to the Beta ones
        for j in range(minBrun,min(i+1,maxBhistory - 1)):
            matrix[i][1][j] = z * matrix[i-1][1][j-1] + x * prxati[i][1] + y * prbrun[j] - y * prbrun[j-1] + w * transition[4][j]
        if i >= maxBhistory - 1:
            matrix[i][1][maxBhistory - 1] = x * sum([prxati[j][1] for j in range(i-maxBhistory+minBrun,i)]) + max([matrix[i - maxBhistory+minBrun][1][l] * z + y * prbrun[maxBhistory - 1] - y * prbrun[l] for l in range(minBrun - 1,maxBhistory)]) + w * sum([transition[4][m] for m in range(minBrun,maxBhistory)])
        if i >= minBrun:
            matrix[i][1][minBrun - 1] = x * sum([prxati[l][1] for l in range(i-minBrun + 1,i + 1)]) + z * max(max([matrix[i-minBrun][0][l] + w * transition[1][l] for l in range(minArun-1,maxAhistory)]),max([matrix[i-minBrun][2][m] + w * transition[7][m] for m in range(minCrun-1,maxChistory)])) + y * prbrun[minBrun - 1] + w * sum([transition[4][m] for m in range(minBrun - 1)])
        if i == minBrun - 1:
            matrix[i][1][minBrun-1] = x * sum([prxati[j][1] for j in range(minBrun)]) + y * prbrun[minBrun-1]
        
        #now for the Coil ones
        for j in range(minCrun,min(i+1,maxChistory)):
            matrix[i][2][j] = z * matrix[i-1][2][j-1] + x * prxati[i][2] + y * prcrun[j] - y * prcrun[j-1] + w * transition[8][j]
        if i >= maxChistory:
            matrix[i][2][maxChistory - 1] = x * sum([prxati[j][2] for j in range(i-maxChistory+minCrun,i-1)]) + max([matrix[i - maxChistory+minCrun][2][l] * z + y * prcrun[maxChistory - 1] - y * prcrun[l] for l in range(minCrun - 1,maxChistory)]) + w * sum([transition[8][m] for m in range(minCrun,maxChistory)])
        if i >= minCrun:
            matrix[i][2][minCrun - 1] = x * sum([prxati[l][2] for l in range(i-minCrun + 1,i + 1)]) + z * max(max([matrix[i-minCrun][0][l] + w * transition[2][l] for l in range(minArun-1,maxAhistory)]),max([matrix[i-minCrun][1][m] + w * transition[5][m] for m in range(minBrun-1,maxBhistory)])) + y * prcrun[minCrun - 1] + w * sum([transition[8][m] for m in range(minCrun-1 )])
        if i == minCrun - 1:
            matrix[i][2][minCrun-1] = x * sum([prxati[j][2] for j in range(minCrun)]) + y * prcrun[minCrun-1]

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
    w = .1
    x = 1
    y = .1
    z = 1 
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
    prruns = [[log(thing) for thing in prrun] for prrun in prruns]
    maxhistories = [len(prruns[x]) for x in range(len(prruns))]
    neginf = -1000000000

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

def backtrack3(matrix):
    numclasses = len(matrix[0])
    maxhistories = [len(matrix[0][k]) for k in range(numclasses)]
    minrunlens = [3,1,1]
    letteranswer = [' ' for i in range(len(matrix))]
    numberanswer = [0 for i in range(len(matrix))]
    lasts = [max(matrix[-1][k][:]) for k in range(numclasses)]
    letters = ['A','B','C']
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
            for j in range(len(letters)):
                if j == letternum:
                    othervals.append(neginf)
                else:
                    othervals.append(max(matrix[i][j][:]))
            maxothervalue = max(othervals)
            letteranswer[i] = letters[othervals.index(maxothervalue)]
            numberanswer[i] = matrix[i][othervals.index(maxothervalue)][:].index(maxothervalue)

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
            for j in range(len(letters)):
                if j == letternum:
                    othervals.append(neginf)
                else:
                    othervals.append(max(matrix[i][j][:]))
            maxothervalue = max(othervals)
            letteranswer[i] = letters[othervals.index(maxothervalue)]
            numberanswer[i] = matrix[i][othervals.index(maxothervalue)][:].index(maxothervalue)

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
    prxati = getprxati(rand)
    if sys.argv[3] == '8':
        matrix = fillmatrixlength8(transitions)
        answer = backtrack8(matrix)
    else:
        matrix = fillmatrixlength(transitions)
        answer = backtrack3(matrix)
except:
    print(traceback.format_exc())
    a,b,c = sys.exc_info()
    answer = ''

out.write(answer + '\n')
overallanswer.append(answer)


