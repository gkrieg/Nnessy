#!/usr/bin/python

import sys

usage = """Usage: %s file1 file2 file3

file1 is a file whose number of points will be counted and written at the top of the file.
file2 is a file whose number of points will be counted and written at the top of the file.
file3 is a file whose number of points will be counted and written at the top of the file.""" % (sys.argv[0])

if len(sys.argv) < 4:
    print(usage, file=sys.stderr)
    exit()

# process the files
for fileName in sys.argv[1:]:

    oldFile = open(fileName, "r")

    lines = oldFile.readlines()
    count = len(lines)

    oldFile.close()

    newFile = open(fileName, "w")
    newFile.write("# Points: %d\n" % (count))

    for line in lines:
        newFile.write(line)

    del lines[:]
    newFile.close()
