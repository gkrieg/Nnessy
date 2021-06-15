import sys

ain = open('23alpha.words','r').readlines()
bein = open('23beta.words','r').readlines()
cin = open('23other.words','r').readlines()
ass = open('alphaSecondaryStructure.words','r').readlines()
bss = open('betaSecondaryStructure.words','r').readlines()
css = open('otherSecondaryStructure.words','r').readlines()
adict = {}
bdict = {}
cdict = {}
aout = open('alpha.words','w')
bout = open('beta.words','w')
cout = open('other.words','w')

def fixss(word1,word2):
    letters = []
    for i in range(len(word1)):
        if word1[i] == word2[i]:
            letters.append(word1[i])
        else:
            if word1[i] == 'A':
                if word2[i] == 'B':
                    letters.append('D')
                elif word2[i] == 'C':
                    letters.append('E')
                elif word2[i] == 'F':
                    letters.append('G')
                else:
                    letters.append(word2[i])
            elif word1[i] == 'B':
                if word2[i] == 'A':
                    letters.append('D')
                elif word2[i] == 'C':
                    letters.append('F')
                elif word2[i] == 'E':
                    letters.append('G')
                else:
                    letters.append(word2[i])
            elif word1[i] == 'C':
                if word2[i] == 'A':
                    letters.append('E')
                elif word2[i] == 'B':
                    letters.append('F')
                elif word2[i] == 'D':
                    letters.append('G')
                else:
                    letters.append(word2[i])
            else:
                letters.append(word1[i])
    return ''.join(letters)



aset = set()
bset = set()
cset = set()

allset = set()

for i, line in enumerate(ain):
    if line not in aset and line not in allset:
        adict[line] = ass[i]
        aset.add(line)
        allset.add(line)
    else:
        adict[line] = fixss(adict[line],ass[i])

for i,line in enumerate(bein):
    if line not in bset and line not in allset:
        bdict[line] = bss[i]
        bset.add(line)
        allset.add(line)
    elif line in bset and line not in aset:
        bdict[line] = fixss(bdict[line],bss[i])
    elif line not in bset and line in aset:
        newword = fixss(adict[line],bss[i])
        adict[line] = newword
        bdict[line] = newword
        bset.add(line)


for i,line in enumerate(cin):
    if line not in cset and line not in allset:
        cdict[line] = css[i]
        cset.add(line)
        allset.add(line)
    elif line in cset and line not in aset and line not in bset:
        cdict[line] = fixss(cdict[line],css[i])
    elif line not in cset and line in allset:
        if line not in bset and line in aset:
            newword = fixss(adict[line],css[i])
            adict[line] = newword
            cdict[line] = newword
        elif line not in aset and line in bset:
            newword = fixss(bdict[line],css[i])
            bdict[line] = newword
            cdict[line] = newword
        elif line in aset and line in bset:
            newword = fixss(bdict[line],css[i])
            adict[line] = newword
            bdict[line] = newword
            cdict[line] = newword
        cset.add(line)

for key in adict:
    aout.write(key.strip() + '_' + adict[key].strip() + '\n')
for key in bdict:
    bout.write(key.strip() + '_' + bdict[key].strip() + '\n')
for key in cdict:
    cout.write(key.strip() + '_' + cdict[key].strip() + '\n')
        
