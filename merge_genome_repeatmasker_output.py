import sys
import argparse
import refMEUtils
# SETUP

parser = argparse.ArgumentParser(description='merge SINE/LINE genome annotations')


parser.add_argument('--repmask', type=str,help='RepeatMasker .out file',required=True)
parser.add_argument('--outsine', type=str,help='output file for merged SINEs',required=True)
parser.add_argument('--outline', type=str,help='output file for merged LINEs',required=True)
args = parser.parse_args()

#####################################################################

print('Processing RepeatMask File',args.repmask)
print('SINE output:',args.outsine)
print('LINE output:',args.outline)


#Read in all RepeatMasker line

allLines = []
inFile = open(args.repmask,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    if len(line) == 0:
        continue
    if line[0] == 'SW':
        continue
    if line[0] == 'score':
        continue
    d = refMEUtils.parse_rmask_line(line)
    
    if d['chrom'] != 'chr1':
        break
    
    allLines.append(d)
inFile.close()
print('Read in all RepeatMasker lines, num lines = ',len(allLines))

print('Merging elements interrupted by Simple_repeats')

toSkip = {}
numMerge = 0
numNotAdj = 0
mergeCounts = {}
for i in range(0,len(allLines)-2):
    if i in toSkip:
        continue
    if allLines[i]['rID'] != allLines[i+2]['rID']:  # want same repeat with something in middle
        continue
    if allLines[i]['ori'] != allLines[i+2]['ori']:  # want same repeat with same orientation
        continue
    if allLines[i]['rFam'] == 'Simple_repeat':   # do not merge interrupted simple repeats
        continue
    if allLines[i+1]['rFam'] != 'Simple_repeat':  # middle must be simple repeat
        continue        
        
     # require that ends be adjacent        
    if (int(allLines[i]['rEnd'])+1) != int(allLines[i+1]['rStart']):
        numNotAdj += 1
        continue

    if (int(allLines[i+2]['rStart'])-1) != int(allLines[i+1]['rEnd']):
        numNotAdj += 1
        continue

    # if get here, can do a merge!
    # change coords of i, add i+1 and i+2 to skip
    numMerge += 1     
    allLines[i]['rEnd'] = allLines[i+2]['rEnd']
     
    allLines[i]['elemStart'] = min(allLines[i]['elemStart'],allLines[i+2]['elemStart'])
    allLines[i]['elemEnd'] = max(allLines[i]['elemEnd'],allLines[i+2]['elemEnd'])     
    allLines[i]['elemLeft'] = max(allLines[i]['elemLeft'],allLines[i+2]['elemLeft'])     
     
    # so we do not consider anymore...
    toSkip[i+1] = 1
    toSkip[i+2] = 1

    if allLines[i]['rFam'] not in mergeCounts:
        mergeCounts[allLines[i]['rFam']] = 0
    mergeCounts[allLines[i]['rFam']] += 1
            
print('\ndid merge',numMerge,len(toSkip))
print('numNotAdj',numNotAdj)
fams = list(mergeCounts.keys())
fams.sort()
for f in fams:
    print(f,mergeCounts[f])

# make new list with just the merges
mergedLines = []
for i in range(0,len(allLines)):
    if i in toSkip:
        continue
    mergedLines.append(allLines[i])
print('merged DONE!',len(mergedLines))

print('\nscanning to merge adjacent LINEs in 2 fragments')

toSkip = {}
numMerge = 0
moreThanTwoFrags = 0
numMergeOriDiff = 0
for i in range(0,len(mergedLines)-1):
    if i in toSkip:
        continue
    if mergedLines[i]['rID'] != mergedLines[i+1]['rID']:  # want same repeat with something in middle      
        continue
    # only do SINE and LINE
    if 'LINE' not in mergedLines[i]['rFam']:
        continue

    numEntries = 0
    searchStart = max(0,i-20)
    searchEnd = min(i+20,len(mergedLines))
    for j in range(searchStart,searchEnd):
        if mergedLines[j]['rID'] == mergedLines[i]['rID']:
            numEntries +=1
    if numEntries != 2:
        moreThanTwoFrags += 1
        continue
    
    if mergedLines[i]['ori'] == mergedLines[i+1]['ori']:  # same orientation, so merge..
    
        numMerge += 1     
        mergedLines[i]['rEnd'] = mergedLines[i+1]['rEnd']     
        mergedLines[i]['elemStart'] = min(mergedLines[i]['elemStart'],mergedLines[i+1]['elemStart'])
        mergedLines[i]['elemEnd'] = max(mergedLines[i]['elemEnd'],mergedLines[i+1]['elemEnd'])     
        mergedLines[i]['elemLeft'] = max(mergedLines[i]['elemLeft'],mergedLines[i+1]['elemLeft'])          
        # so we do not consider anymore...
        toSkip[i+1] = 1
    else:  # ori does not match        
        # take orientation of end closes to 3' end
        if mergedLines[i]['elemEnd'] > mergedLines[i+1]['elemEnd']:
            newOri = mergedLines[i]['ori']
        else:
            newOri = mergedLines[i+1]['ori']
        numMergeOriDiff += 1     
        mergedLines[i]['rEnd'] = mergedLines[i+1]['rEnd']     
        mergedLines[i]['elemStart'] = min(mergedLines[i]['elemStart'],mergedLines[i+1]['elemStart'])
        mergedLines[i]['elemEnd'] = max(mergedLines[i]['elemEnd'],mergedLines[i+1]['elemEnd'])     
        mergedLines[i]['elemLeft'] = max(mergedLines[i]['elemLeft'],mergedLines[i+1]['elemLeft'])   
        mergedLines[i]['ori'] = newOri
        # so we do not consider anymore...
        toSkip[i+1] = 1

print('\ndid merge',numMerge,len(toSkip))
print('numMergeOriDiff',numMergeOriDiff)

# make new list with just the merges
mergedLines2 = []
for i in range(0,len(mergedLines)):
    if i in toSkip:
        continue
    mergedLines2.append(mergedLines[i])
print('merged2 DONE!',len(mergedLines2))

print('ready for output!')

SINE_locus_count = 0
LINE_locus_count = 0

outSINE = open(args.outsine,'w')
outLINE = open(args.outline,'w')

nl = ['#locusID','chrom','rStart','rEnd','ori','rType','rFam','elemStart','elemEnd','elemLeft','dvg','rIDs']
nl = '\t'.join(nl) + '\n'
outSINE.write(nl)
outLINE.write(nl)

for i in range(0,len(mergedLines2)):
    nl = []
    d = mergedLines2[i]
    if 'SINE' in d['rFam']:
        lid = f'SINE_{SINE_locus_count}'
        SINE_locus_count += 1
        nl.append(lid)
        nl.append(d['chrom'])
        nl.append(d['rStart'])
        nl.append(d['rEnd'])
        nl.append(d['ori'])
        nl.append(d['rType'])
        nl.append(d['rFam'])
        nl.append(d['elemStart'])
        nl.append(d['elemEnd'])
        nl.append(d['elemLeft'])
        nl.append(d['dvg'])
        nl.append(d['rID'])        
        nl = [str(j) for j in nl]
        nl = '\t'.join(nl) + '\n'
        outSINE.write(nl)    
    elif 'LINE' in d['rFam']:
        lid = f'LINE_{LINE_locus_count}'
        LINE_locus_count += 1
        nl.append(lid)
        nl.append(d['chrom'])
        nl.append(d['rStart'])
        nl.append(d['rEnd'])
        nl.append(d['ori'])
        nl.append(d['rType'])
        nl.append(d['rFam'])
        nl.append(d['elemStart'])
        nl.append(d['elemEnd'])
        nl.append(d['elemLeft'])
        nl.append(d['dvg'])
        nl.append(d['rID'])        
        
        nl = [str(j) for j in nl]
        nl = '\t'.join(nl) + '\n'
        outLINE.write(nl)    
outSINE.close()
outLINE.close()

print('DONE!')
