import sys
import argparse
import refMEUtils

#####################################################################################
# SETUP
parser = argparse.ArgumentParser(description='filter for LINE/L1 and 3\' end')



parser.add_argument('--elementtable', type=str,help='input file element table',required=True)
parser.add_argument('--parsetable', type=str,help='parsed align table',required=True)
parser.add_argument('--tabulattable', type=str,help='tabulated output table',required=True)

args = parser.parse_args()
#####################################################################################



print('Reading in TSD position table',args.elementtable)

alleleTable = {}
inFile = open(args.elementtable,'r')
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    lid = line[0]
    empty = line[1]
    filled = line[2]
    empty = empty.split('-')
    eS = int(empty[0])
    eE = int(empty[1])
    
    filled = filled.split(',')
    fL = filled[0]
    fL = fL.split('-')
    fLS = int(fL[0])
    fLE = int(fL[1])

    fR = filled[1]
    fR = fR.split('-')
    fRS = int(fR[0])
    fRE = int(fR[1])
    
    alleleTable[lid] = [ [eS,eE], [fLS,fLE],[fRS,fRE]]
inFile.close()

print('read in %i allele info!' % len(alleleTable))

counts = {}
for lid in alleleTable:
    counts[lid] = {}
    counts[lid]['empty'] = 0
    counts[lid]['filled5'] = 0
    counts[lid]['filled3'] = 0    
    
print('setup counts table',len(counts))

inFile = open(args.parsetable,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    lid = line[0]
    allele = line[1]
    res = line[2]
    
    if allele == 'empty' and res == 'PRESENT':
        counts[lid]['empty'] += 1
    elif allele == 'filled':
        res = res.split(',')
        if res[0] == '5prime_PRESENT':
            counts[lid]['filled5'] += 1
        if res[1] == '3prime_PRESENT':
            counts[lid]['filled3'] += 1
inFile.close()

print('tabulated all, writing output!')
outFile = open(args.tabulattable,'w')
for lid in alleleTable:
    nl = [lid,counts[lid]['empty'],counts[lid]['filled5'],counts[lid]['filled3']]
    
    dec = '?'
    sumFilled = counts[lid]['filled5'] + counts[lid]['filled3']
    if counts[lid]['empty'] > 0 and sumFilled > 0:
        dec = 'obs_both'
    elif counts[lid]['empty'] == 0 and sumFilled == 0:
        dec = 'obs_neither'
    elif counts[lid]['empty'] > 1 and sumFilled == 0:
        dec = 'mult_empty'
    elif counts[lid]['empty'] == 1 and sumFilled == 0:
        dec = 'EMPTY'
    elif counts[lid]['filled5'] > 1 and counts[lid]['empty'] == 0:
        dec = 'mult_filled'
    elif counts[lid]['filled3'] > 1 and counts[lid]['empty'] == 0:
        dec = 'mult_filled'
    elif counts[lid]['filled5'] == 1 and counts[lid]['filled3'] == 0 and counts[lid]['empty'] == 0:
        dec = 'FILLED_5'
    elif counts[lid]['filled3'] == 1 and counts[lid]['filled5'] == 0 and counts[lid]['empty'] == 0:
        dec = 'FILLED_3'
    elif counts[lid]['filled5'] == 1 and counts[lid]['filled3'] == 1 and counts[lid]['empty'] == 0:
        dec = 'FILLED_both'
    else:
        print(nl)
        print('what??')
        sys.exit()

    nl.append(dec)     

    
    
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
outFile.close()
print('done!')
