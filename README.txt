To install Nnessy, run:

./install


To run Nnessy, make sure that your input file is in .fasta format and contains a single sequence


./nnessy -i [input].fasta [-n numthreads] [-o output]


If no output file is specified, the secondary structure will be placed in a file with the name output/[input].ss, which has the secondary structure.

-n option allows the number of threads to be specified if run in a parallel setting.

