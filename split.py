import sys

wordLength = 23
halfLength = wordLength//2
def sortSequence(sequence):
    allwords = []
    for i in range(len(sequence)):
        if i < wordLength//2 :
            word = ''.join(['Z' for j in range( halfLength - i)]) + sequence[0:(i + halfLength + 1)]
        elif i >= len(sequence) - 1 - halfLength:
            fromend = len(sequence) - i -1
            word = sequence[i - halfLength:] + ''.join(['Z' for j in range( halfLength - fromend)])
        else:
            word = sequence[i - halfLength : i + halfLength + 1]
        containsLower = False
        for c in word:
            if c.islower():
                containsLower = True
        if not containsLower:
            allwords.append(word)
    return allwords

infile = open(sys.argv[1],'r')

inlines = infile.readlines()
numseqs = len(inlines)  

for i in range(1,numseqs,2):
    outfile = open(sys.argv[1] + '.seq','w')
    allwords = sortSequence(inlines[i].strip())
    outfile.write('# Points: {0}\n'.format(len(allwords)))
    for word in allwords:
        outfile.write(word + '\n')

