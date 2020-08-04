import sys

def func(x,a,b):
    return a * x + b

def getdistacc(a,b,c,d,val,e):
    return func(val*100,a,b)/100 if val * 100 < e else func(val*100,c,d)/100

def getestimatedist(f,prefix,returnarray=False):
    if prefix == '8state':
        a = -0.19800842564697718 
        b = 98.28064203538942 
        c = -6.309037446564571 
        d = 447.1813308129393
        e = 57
    else:
        a = -0.824600555782599 
        b = 133.42144454707145 
        c = -5.855914499014564 
        d = 432.5124697999706
        e = 61
    accs = []
    for val in f:
        val = float(val.strip())
        boxacc = getdistacc(a,b,c,d,val,e)
        accs.append(boxacc)
    if returnarray == True:
        return accs
    else:
        return sum(accs)/len(accs)

accs = getestimatedist(open(sys.argv[1],'r').readlines(),sys.argv[2])

output = open(sys.argv[3],'w')
output.write(str(accs))
