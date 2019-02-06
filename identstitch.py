import sys
iin = open(sys.argv[1],'r').readlines()
sin = open(sys.argv[2],'r').readlines()
outf = open(sys.argv[3],'w')
k = int(sys.argv[4])

i = 0
iinpointer = 0
sinpointer = 0
endprotein = ((len(iin) + len(sin)) // k ) - 1
while i < endprotein:
    print('len',((len(iin) + len(sin)) // k) - 1,i)
    ii = int(iin[iinpointer ].split()[0]) if iinpointer < len(iin) else endprotein
    print(i,ii)
    while i < ii:
        print(i)
        for l in range(k):
            print(sin[sinpointer])
            outf.write(sin[sinpointer])
            sinpointer += 1
        i += 1
    if ii != endprotein:
        for l in range(k):
            print(iin[iinpointer])
            outf.write(iin[iinpointer])
            iinpointer += 1
    i += 1
