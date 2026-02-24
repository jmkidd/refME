import sys
import argparse
import refMEUtils
from Bio import Align

# SETUP

parser = argparse.ArgumentParser(description='annotate alignments for TSD and poly(A)')


parser.add_argument('--fasta', type=str,help='genome fata file',required=True)
parser.add_argument('--intable', type=str,help='input file of segments',required=True)
parser.add_argument('--outtable', type=str,help='output file for segments with annotation',required=True)
args = parser.parse_args()




#####################################################################
def search_tsd_polyA(d,chromSeq):
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match = 1
    aligner.mismatch = -2
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    
    chromLen = len(chromSeq)
    
    five_delta = 100
    five_internal = 50
    three_internal = 50
    three_delta = 3000
    
    if d['ori'] == '+':
        fiveStart = d['rStart'] - five_delta
        fiveEnd = d['rStart'] + five_internal
        
        threeStart = d['rEnd'] - three_internal
        threeEnd = d['rEnd'] + three_delta
        
    
    elif d['ori'] == 'C':
        fiveStart = d['rEnd'] - five_internal
        fiveEnd = d['rEnd'] + five_delta
        
        threeStart = d['rStart'] - three_delta
        threeEnd = d['rStart'] + three_internal
    else:
        print('orient error?')
        print(d)
        sys.exit()
        
        
    print(fiveStart,fiveEnd,threeStart,threeEnd)
    
    if fiveStart < 1:
        fiveStart = 1
    if threeStart < 1:
        threeStart = 1
    
    if fiveEnd > chromLen:
        fiveEnd = chromLen
    if threeEnd > chromLen:
        threeEnd = chromLen
    
    
    fiveSeq = chromSeq[fiveStart-1:fiveEnd]
    threeSeq = chromSeq[threeStart-1:threeEnd]
    
    alignments = aligner.align(fiveSeq, threeSeq)

    print('num alignments',len(alignments))

    annotCandidates = []
    for a_i, alignment in enumerate(alignments):
        print(a_i)
        print(alignment)
        
        
        aln_five = alignment[0]
        aln_three = alignment[1]
        
        aln_five_coord = alignment.coordinates[0]
        aln_three_coord = alignment.coordinates[1]        

        aln_five_start = aln_five_coord[0] + 1
        aln_five_end = aln_five_coord[1]
        
        aln_three_start = aln_three_coord[0] + 1
        aln_three_end = aln_three_coord[1]
        
        
        aln_five_start = aln_five_start + fiveStart - 1
        aln_five_end = aln_five_end + fiveStart - 1
        

        aln_three_start = aln_three_start + threeStart - 1
        aln_three_end = aln_three_end + threeStart - 1
        
        
        
        
        
        
        print(aln_five,aln_three)
        print(aln_five_coord,aln_three_coord)
        print(aln_five_start,aln_five_end,aln_three_start,aln_three_end)
        candidate = score_candidate(aln_five,aln_three,aln_five_start,aln_five_end,aln_three_start,aln_three_end,d,chromSeq)
        annotCandidates.append(candidate)
                
    print('there are ',len(annotCandidates))
    maxScore = -1 
    for i in range(len(annotCandidates)):
        if annotCandidates[i]['TSDscore'] > maxScore:
            maxScore = annotCandidates[i]['TSDscore']
        
    if maxScore <= 0: # no alignment, or none with TSD that passess
        res = {}
        res['aln_five'] = '.'
        res['aln_three'] = '.'
        res['aln_five_start'] = '.'
        res['aln_five_end'] = '.'
        res['aln_three_start'] = '.'
        res['aln_three_end'] = '.'        
        res['tsdLen'] = 0
        res['TSDscore'] = 0
        res['tsd_dist_5'] = '.'
        res['tsd_dist_3'] = '.'
        res['passWinPolyA'] = 'NO'
        return res
    
    canAligns = []
    for i in range(len(annotCandidates)): 
        if annotCandidates[i]['TSDscore'] == maxScore:
            canAligns.append(annotCandidates[i])
    
    print('max score is',maxScore,'num with max is',len(canAligns))
    # see how many
    if len(canAligns) == 1:
        return canAligns[0]
    # there are multiple with same score...
    # prefer passWinPolyA yes, after that take smallest sum of distance....
    canAligns = []
    for i in range(len(annotCandidates)): 
        if annotCandidates[i]['TSDscore'] == maxScore and annotCandidates[i]['passWinPolyA'] == 'YES':
            canAligns.append(annotCandidates[i])
    
    print('num with passPolyA is',len(canAligns))
    
    if len(canAligns) == 1:
        return canAligns[0]
    if len(canAligns) == 0: # none of the top have polyA, so need to redo, if there are mult, will take...        
        canAligns = []
        for i in range(len(annotCandidates)): 
            if annotCandidates[i]['TSDscore'] == maxScore:
                canAligns.append(annotCandidates[i])
    print('pick based on dist',len(canAligns))

    #will now take ones with min dist, or just first one
    minDist = 9999999999
    for i in range(len(canAligns)): 
        dist = canAligns[i]['tsd_dist_5'] + canAligns[i]['tsd_dist_3']
        if dist < minDist:
            minDist = dist
    print('mindist is',minDist)
    # return the min
    for i in range(len(canAligns)): 
        dist = canAligns[i]['tsd_dist_5'] + canAligns[i]['tsd_dist_3']
        if dist == minDist:
            return(canAligns[i])    
#####################################################################
def score_candidate(aln_five,aln_three,aln_five_start,aln_five_end,aln_three_start,aln_three_end,d,chromSeq):
    print('in score!!')
    
    myCan = {} # store candidate data
    myCan['aln_five'] = aln_five
    myCan['aln_three'] = aln_three
    myCan['aln_five_start'] = aln_five_start
    myCan['aln_five_end'] = aln_five_end
    myCan['aln_three_start'] = aln_three_start
    myCan['aln_three_end'] = aln_three_end
    
    fiveLen = myCan['aln_five_end'] - myCan['aln_five_start'] + 1
    threeLen = myCan['aln_three_end'] - myCan['aln_three_start'] + 1
    
    myCan['tsdLen'] = min(fiveLen,threeLen)
    
    
    numMatch = 0
    numMissMatch = 0
    for i in range(len(myCan['aln_five'])):
        if myCan['aln_five'][i] == myCan['aln_three'][i]:
            numMatch += 1
        else:
            numMissMatch += 1
            
    myCan['numMatch'] = numMatch
    myCan['numMissMatch'] = numMissMatch
    
    
    print(fiveLen,threeLen)
    print('matches',numMatch,numMissMatch)
    
    if '-' in myCan['aln_five'] or '-' in myCan['aln_three']:
        print('true gap!')
        myCan['hasGap'] = 'YES'
    else:
        print('no gap')
        myCan['hasGap'] = 'NO'
    
    
    # do some scoring
    if myCan['tsdLen'] < 9:
        myCan['TSDscore'] = 0
        myCan['tsd_dist_5'] = '.'
        myCan['tsd_dist_3'] = '.'
        myCan['passWinPolyA'] = 'NO'
        myCan['cutSite'] = '.'
        return myCan
        
    # follow rule mosty based on  https://pubmed.ncbi.nlm.nih.gov/12372140/   
    HSP_position_score = 100
    
    if d['ori'] == '+':
        elemStart = d['rStart']
        TSDedge = myCan['aln_five_end']
        tsdDist = d['rStart']- myCan['aln_five_end'] 
    else:    
        elemStart = d['rEnd']
        TSDedge = myCan['aln_five_start']
        tsdDist = TSDedge - elemStart
    print('5tsd dist',tsdDist)
    myCan['tsd_dist_5'] = tsdDist
    
    if tsdDist <= 10:
        HSP_position_score_5 = 100
    else:
        HSP_position_score_5 = 100 - tsdDist*(8.0/9.0)
        
        
    if d['ori'] == '+':
        elemEnd = d['rEnd']
        TSDedge = myCan['aln_three_start']
        tsdDist = TSDedge - elemEnd
    else:    
        elemStart = d['rStart']
        TSDedge = myCan['aln_three_end']
        tsdDist = elemStart - TSDedge
    print(TSDedge,'3tsd dist',tsdDist)
    myCan['tsd_dist_3'] = tsdDist
    
    
    if tsdDist <= 30:
        HSP_position_score_3 = 20
    else:
        HSP_position_score_3 = 0
        
    if myCan['tsdLen'] >= 11 and myCan['tsdLen'] <= 18 and myCan['numMissMatch'] == 0:
        HSP_qual_score = 20
    elif myCan['tsdLen'] < 11 and myCan['tsdLen'] >= 9  and myCan['numMissMatch'] == 0:
        HSP_qual_score = 10
    elif myCan['tsdLen'] <= 18 and myCan['numMissMatch'] <= 1:
        HSP_qual_score = 10
    elif myCan['tsdLen'] > 18 and myCan['numMissMatch'] <= (myCan['tsdLen'] /50):
        HSP_qual_score = 10
    else:
        HSP_qual_score = 0
        
    myCan['TSDscore'] = HSP_position_score_5 + HSP_position_score_3 + HSP_qual_score
    print('SCORE:',myCan['TSDscore'],HSP_position_score_5,HSP_position_score_3,HSP_qual_score)
    
    print('SEARCH poly(A)!')
    # min 10, at least 73% A...
    # matt took 30bp window, required it to have at least 1 window with 10/15 A
    myCan['passWinPolyA'] = 'NO'
    if d['ori'] == '+':
        endBeforeTSD = myCan['aln_three_start'] -1
        seqStart = endBeforeTSD - 30 + 1
        winSeq = chromSeq[seqStart-1:endBeforeTSD]
        winSeq = winSeq.upper()
        for i in range(0,30-15):
            s = winSeq[i:i+15]
            numA = s.count('A')
            if numA >= 10:
                myCan['passWinPolyA'] = 'YES'
        
    else:
        endBeforeTSD = myCan['aln_three_end'] +1
        seqStart = endBeforeTSD + 30 - 1
        winSeq = chromSeq[endBeforeTSD-1:seqStart]
        winSeq = winSeq.upper()
        for i in range(0,30-15):
            s = winSeq[i:i+15]
            numA = s.count('T')
            if numA >= 10:
                myCan['passWinPolyA'] = 'YES'
        
    
    # get EN cutsite, 5+2
    
    if d['ori'] == '+':
        cutStart = myCan['aln_five_start'] - 2
        cutEnd = myCan['aln_five_start'] + 4
        cutSeq = chromSeq[cutStart-1:cutEnd]
        cutSeq = refMEUtils.revcomp(cutSeq)    
        myCan['cutSite'] = cutSeq
    else:
        cutEnd = myCan['aln_five_end'] + 2
        cutStart = myCan['aln_five_end'] - 4
        cutSeq = chromSeq[cutStart-1:cutEnd]
        myCan['cutSite'] = cutSeq
    
    
    
    return(myCan)

#####################################################################
print('input table:',args.intable)
print('output table:',args.outtable)

# read in sequences

genomeSeqs = refMEUtils.read_fasta_file_to_dict(args.fasta)
n = list(genomeSeqs.keys())
numSeqs = len(n)
print(f'read in {numSeqs} from fasta file')


numLoci = 0
outFile = open(args.outtable,'w')
inFile = open(args.intable,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    if line[0][0] == '#': # header
        print('header!')
        print(line)
        headerRow = []
        colToNames = {}
        for i,n in enumerate(line):
            if n[0] == '#':
                n=n[1:]                
            colToNames[n] = i
            headerRow.append(n)
        
        nl = []
        for j in headerRow:
            nl.append(j)
        nl[0] = '#' + nl[0]
        
        nl.append('tsdLen')
        nl.append('TSDscore')
        nl.append('passWindowPolyA')
        nl.append('TSD-5prime-coords')
        nl.append('TSD-5prime-seq')
        nl.append('TSD-3prime-coords')
        nl.append('TSD-3prime-seq')
        nl.append('TSD-5prime-distance')
        nl.append('TSD-3prime-distance')        
        nl.append('inferredCleavageSeq')
        


        
        nl = '\t'.join(nl) + '\n'
        outFile.write(nl)
        
        continue
    # not header
    print(line)
    d = {}
    for n in colToNames:
        v = line[colToNames[n]]
        if n in ['rStart','rEnd','elemStart','elemEnd']:
            v = int(v)
        d[n] = v
    print(d)
    
    chromSeq = genomeSeqs[d['chrom']]['seq']
    
    res = search_tsd_polyA(d,chromSeq)
    print('my result is!!!')
    print(res)
    
    
    nl = []
    for j in headerRow:
        nl.append(d[j])
        
    
    nl.append(res['tsdLen'])    
    nl.append(res['TSDscore'])    
    nl.append(res['passWinPolyA'])
    
    
    fiveTSDCoords = '%i-%i' % (res['aln_five_start'],res['aln_five_end'])
    threeTSDCoords = '%i-%i' % (res['aln_three_start'],res['aln_three_end'])
    
    
    nl.append(fiveTSDCoords)
    
    tsd  = res['aln_five'].replace('-','')
    if d['ori'] == 'C':
        tsd = refMEUtils.revcomp(tsd)
    nl.append(tsd)
    
    nl.append(threeTSDCoords)

    tsd  = res['aln_three'].replace('-','')
    if d['ori'] == 'C':
        tsd = refMEUtils.revcomp(tsd)
    nl.append(tsd)


    nl.append(res['tsd_dist_5'])
    nl.append(res['tsd_dist_3'])
    nl.append(res['cutSite'])
        
    
    nl = [str(j) for j in nl]    
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
    
    numLoci += 1
    if numLoci >= 100:
        break
        
    
    



inFile.close()
outFile.close()
