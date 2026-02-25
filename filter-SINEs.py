import sys
import argparse

# SETUP

parser = argparse.ArgumentParser(description='filter for LINE/L1 and 3\' end')



parser.add_argument('--intable', type=str,help='input file',required=True)
parser.add_argument('--outtable', type=str,help='output file',required=True)

args = parser.parse_args()


inFile = open(args.intable,'r')
outFile = open(args.outtable,'w')
for line in inFile:
    ol = line
    line = line.rstrip()
    line = line.split()
    if line[0][0] == '#': # header
        outFile.write(ol)
        print('header!')
        print(line)
        headerRow = []
        colToNames = {}
        for i,n in enumerate(line):
            if n[0] == '#':
                n=n[1:]                
            colToNames[n] = i
            headerRow.append(n)
        continue
    if line[colToNames['rFam']] != 'SINE/tRNA':
        continue
   
    rType = line[colToNames['rType']]

    if rType in ['LFSINE_Vert']:
        continue    
   
    
    
    
    outFile.write(ol)
    
inFile.close()
outFile.close()

