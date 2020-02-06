
import sys
rand = sys.argv[1]
outfile = open('allgreedy{}'.format(rand),'w')

for i in range(int(sys.argv[2])):
    try:
        infile = open('predictions/{0}/probabilities{1}'.format(rand,i),'r').readlines()
        #infile = open('predictions/{0}/newprobs{1}'.format(rand,i),'r').readlines()
        #infile = open('predictions/{0}/cprobs{1}'.format(rand,i),'r').readlines()
        #infile = open('outtt','r').readlines()

        for line in infile:
            a = float(line.split()[0].strip())
            b = float(line.split()[1].strip())
            c = float(line.split()[2].strip())
            if a > b and a > c:
                ans = 'A'
            elif b > c:
                ans = 'B'
            else:
                ans = 'C'
            outfile.write('{0}'.format(ans))
    except:
        print('oops {0}'.format(i))
    outfile.write('\n')
