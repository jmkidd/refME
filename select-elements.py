import sys
import argparse

# SETUP

parser = argparse.ArgumentParser(description='selected elements based on 5\' TSD location, also filter N and chrom')



parser.add_argument('--intable', type=str,help='input file',required=True)
parser.add_argument('--outtable', type=str,help='output file',required=True)

args = parser.parse_args()

MIN_TSD_5_DIST = -10
MAX_TSD_5_DIST = 20


chromToDo = [f'chr{c}' for c in range(1,39)]
chromToDo.append('chrX')

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

    if 'N' in line[colToNames['TSD-5prime-seq']]:
        continue
    if 'N' in line[colToNames['TSD-3prime-seq']]:
        continue
    if 'N' in line[colToNames['inferredCleavageSeq']]:
        continue
        
    tsdLen = int(line[colToNames['tsdLen']])
    if tsdLen < 9:
        continue
        
    if line[colToNames['passWindowPolyA']] == 'NO':
        continue
    
    
        
    tsd5dist = int(line[colToNames['TSD-5prime-distance']])
    
    if tsd5dist < MIN_TSD_5_DIST or tsd5dist > MAX_TSD_5_DIST:
        continue

    if line[colToNames['chrom']] not in chromToDo:
        continue


    outFile.write(ol)
    
inFile.close()
outFile.close()

