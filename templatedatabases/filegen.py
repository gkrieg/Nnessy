import sys
from array import array

def simplifyAlphabet(word):
    carray = array('b', word.encode())
    for c in range(len(word)):
        if word[c] == 'H' or word[c] == 'h' or word[c] == 'g' or word[c] == 'G' or word[c] == 'I' or word[c] == 'i':
            carray[c] = ord('A')
        elif word[c] == 'E' or word[c] == 'e' or word[c] == 'B' or word[c] == 'b':
            carray[c] = ord('B')
        else:
            carray[c] = ord('C')
    word = carray.tobytes().decode()
    return word


def addZs(word,halflen):
    arr = ['Z' for x in range(halflen)] + [y for y in word] + ['Z' for z in range(halflen)]
    return ''.join(arr)

def refineaastructure(word):

    residues = "ARNDCQEGHILKMFPSTWYVBX"
    carray = array('b', word.encode())
    for c in range(len(word)):
        if word[c] not in residues:
            carray[c] = ord('X')
        else:
            carray[c] = ord(word[c])
    word = carray.tobytes().decode()
    return word


inf = open(sys.argv[1],'r').readlines()

wordlen = int(sys.argv[2])
aout = open('{}alpha.words'.format(wordlen),'w')
asout = open('alphaSecondaryStructure.words','w')
bout = open('{}beta.words'.format(wordlen),'w')
bsout = open('betaSecondaryStructure.words','w')
cout = open('{}other.words'.format(wordlen),'w')
csout = open('otherSecondaryStructure.words','w')

halflen = wordlen // 2

for i in range(0,len(inf),5):
    relnum = i // 5
    aa = inf[i+1].strip()
    ss = inf[i+3].strip()
    ss = simplifyAlphabet(ss)
    aa = refineaastructure(aa)
    aa = addZs(aa,halflen)
    ss = addZs(ss,halflen)

    for j in range(len(aa) - wordlen):
        aaword = aa[j:j+wordlen]
        ssword = ss[j:j+wordlen]
        if ssword[halflen] == 'A':
            aout.write(aaword + '\n')
            asout.write(ssword + '#{}'.format(relnum)+ '\n')
        elif ssword[halflen] == 'B':
            bout.write(aaword + '\n')
            bsout.write(ssword + '#{}'.format(relnum)+ '\n')
        elif ssword[halflen] == 'C':
            cout.write(aaword + '\n')
            csout.write(ssword + '#{}'.format(relnum)+ '\n')


