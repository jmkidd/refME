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
        
        break
        
    print('there are ',len(annotCandidates))
    res = {}    
    return res
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
    
    numMatch = 0
    numMissMatch = 0
    for i in range(len(myCan['aln_five'])):
        if myCan['aln_five'][i] == myCan['aln_three'][i]:
            numMatch += 1
        else:
            numMissMatch += 1
            
    
    
    print(fiveLen,threeLen)
    print('matches',numMatch,numMissMatch)
    
    if '-' in myCan['aln_five'] or '-' in myCan['aln_three']:
        print('true gap!')
    else:
        print('no gap')
    
    
    
    
    
    return(myCan)
    
    
    



#####################################################################
print('input table:',args.intable)
print('output table:',args.outtable)

# read in sequences

genomeSeqs = refMEUtils.read_fasta_file_to_dict(args.fasta)
n = list(genomeSeqs.keys())
numSeqs = len(n)
print(f'read in {numSeqs} from fasta file')


outFile = open(args.outtable,'w')
inFile = open(args.intable,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    if line[0][0] == '#': # header
        print('header!')
        print(line)
        colToNames = {}
        for i,n in enumerate(line):
            if n[0] == '#':
                n=n[1:]
            colToNames[n] = i
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
    
    

    
    
    
    sys.exit()
    



inFile.close()
outFile.close()
