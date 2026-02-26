import sys
import argparse
import refMEUtils

#####################################################################################
def score_aln(paf,cigarExp,alleleTable):
    EMPTY_FLANK_SIZE = 20
    EMPTY_MISS_PERMIT = 2
    
    TSD_DEL_SEARCH_SIZE = 10
    TSD_DEL_PERMIT = 10
    
    
    lid = paf['qName']
    qLen = paf['qLen']
    p = lid.split('_')
    locusName = p[0] + '_' + p[1]
    alleleType = p[2]
    res = {}
    res['locusName'] = locusName
    res['alleleType'] = alleleType
    
    # keep track of which bp in the query are aligned (part of M)
    # do not care of sequence matches or does not
    qMatches = {}
    qAttDels = {}
    for i in range(1,paf['qLen']+1):
        qAttDels[i] = 0
    
    if paf['strand'] == '+':
        currentPos = paf['qStart'] + 1
        for op in cigarExp:
            if op[1] == 'M':
                for i in range(op[0]):
                    qMatches[currentPos] =1
                    currentPos += 1
            elif op[1] == 'I':
                for i in range(op[0]):
                    currentPos += 1
            elif op[1] == 'D':
                qAttDels[currentPos] += op[0]
                    
    else: # on - strand
        currentPos = paf['qEnd']
        for op in cigarExp:
            if op[1] == 'M':
                for i in range(op[0]):
                    qMatches[currentPos] =1
                    currentPos -= 1
            elif op[1] == 'I':
                for i in range(op[0]):
                    currentPos -= 1
            elif op[1] == 'D':
                qAttDels[currentPos] += op[0]
                    
    if res['alleleType'] == 'empty':
        numPresent = 0
        numMiss = 0
        TSDstart = alleleTable[locusName][0][0]
        TSDend = alleleTable[locusName][0][1]
        regStart = TSDstart -1 - EMPTY_FLANK_SIZE + 1
        regEnd = TSDend + 1 + EMPTY_FLANK_SIZE -1
        for i in range(regStart,regEnd+1):
            if i in qMatches:
                numPresent += 1
            else:
                numMiss += 1
        
        numDelRegion = 0
        searchStart = TSDstart -1 -TSD_DEL_SEARCH_SIZE + 1
        searchEnd = TSDend + 1 + TSD_DEL_SEARCH_SIZE -1
        for i in range(searchStart,searchEnd+1):
            numDelRegion += qAttDels[i]
        
        if numDelRegion > TSD_DEL_PERMIT:
            res['decision'] = 'FAIL_del'
        else:        
            if numMiss <= EMPTY_MISS_PERMIT:
                res['decision'] = 'PRESENT'
            else:
                res['decision'] = 'FAIL'
    else: # is filled site 
        numPresentFive = 0
        numMissFive = 0
        TSDstart = alleleTable[locusName][1][0]
        TSDend = alleleTable[locusName][1][1]
        regStart = TSDstart -1 - EMPTY_FLANK_SIZE + 1
        regEnd = TSDend + 1 + EMPTY_FLANK_SIZE -1
        for i in range(regStart,regEnd+1):
            if i in qMatches:
                numPresentFive += 1
            else:
                numMissFive += 1
        
        numDelRegion = 0
        searchStart = TSDstart -1 -TSD_DEL_SEARCH_SIZE + 1
        searchEnd = TSDend + 1 + TSD_DEL_SEARCH_SIZE -1
        for i in range(searchStart,searchEnd+1):
            numDelRegion += qAttDels[i]

        if numDelRegion > TSD_DEL_PERMIT:
            res['decision'] = '5prime_FAIL_del'
        else:
            if numMissFive <= EMPTY_MISS_PERMIT:
                res['decision'] = '5prime_PRESENT'
            else:
                res['decision'] = '5prime_FAIL'
        
        numPresentThree = 0
        numMissThree = 0
        TSDstart = alleleTable[locusName][2][0]
        TSDend = alleleTable[locusName][2][1]
        regStart = TSDstart -1 - EMPTY_FLANK_SIZE + 1
        regEnd = TSDend + 1 + EMPTY_FLANK_SIZE -1
        for i in range(regStart,regEnd+1):
            if i in qMatches:
                numPresentThree += 1
            else:
                numMissThree += 1
        numDelRegion = 0
        searchStart = TSDstart -1 -TSD_DEL_SEARCH_SIZE + 1
        searchEnd = TSDend + 1 + TSD_DEL_SEARCH_SIZE -1
        for i in range(searchStart,searchEnd+1):
            numDelRegion += qAttDels[i]

        if numDelRegion > TSD_DEL_PERMIT:
            res['decision'] += ',3prime_FAIL_del'
        else:
            if numMissThree <= EMPTY_MISS_PERMIT:
                res['decision']+=',3prime_PRESENT'
            else:
                res['decision']+= ',3prime_FAIL'
            
    return res
#####################################################################################

# SETUP
parser = argparse.ArgumentParser(description='filter for LINE/L1 and 3\' end')



parser.add_argument('--elementtable', type=str,help='input file element table',required=True)
parser.add_argument('--paf', type=str,help='paf aligment file',required=True)
parser.add_argument('--outtable', type=str,help='output table',required=True)

args = parser.parse_args()

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


inFile = open(args.paf,'r')
outFile = open(args.outtable,'w')
numDid = 0
for line in inFile:
    line = line.rstrip()
    line = line.split()
    paf = refMEUtils.parse_paf_line(line)
    # get cigar
    cg = '?'
    for i in range(12,len(line)):
        if line[i][0:3]== 'cg:':
            cg = line[i]
            break
    if cg == '?':
        print('no ciga??')
        print(line)
        sys.exit()
    cigar = cg.split(':')[2]
    cigarExp = refMEUtils.expand_cigar(cigar)


    res = score_aln(paf,cigarExp,alleleTable)
    
    nl = []
    nl.append(res['locusName'])
    nl.append(res['alleleType'])
    nl.append(res['decision'])
    nl.append(paf['mapQ'])
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
    
    numDid += 1
    if numDid > 500:    
        break




inFile.close()
outFile.close()