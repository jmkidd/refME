import sys
import argparse
import refMEUtils
#####################################################################
def extract_element(d,chromSeq):
    res = {}
    if d['ori'] == '+':
        res['start'] = d['TSD-5prime-coords'][1] + 1
        res['end'] = d['TSD-3prime-coords'][0] - 1
        res['seq'] = chromSeq[res['start']-1:res['end']]
    else:
        res['start'] = d['TSD-3prime-coords'][1] + 1
        res['end'] = d['TSD-5prime-coords'][0] - 1
        res['seq'] = refMEUtils.revcomp(chromSeq[res['start']-1:res['end']])
    
    res['seq'] = res['seq'].upper()    
    return res
#####################################################################
def get_filled_empty(d,chromSeq,flankSize):
    tsd5Start = d['TSD-5prime-coords'][0]
    tsd5End = d['TSD-5prime-coords'][1]    
    tsd3Start = d['TSD-3prime-coords'][0]
    tsd3End = d['TSD-3prime-coords'][1]    
    
    res = {}
        
    if d['ori'] == '+':
         leftEnd = tsd5Start - 1
         leftStart = leftEnd - flankSize + 1
         midStart = tsd5End + 1
         midEnd = tsd3Start - 1
         rightStart = tsd3End + 1
         rightEnd = rightStart + flankSize - 1
         
         leftSeq = chromSeq[leftStart-1:leftEnd]
         leftTSDSeq = chromSeq[tsd5Start-1:tsd5End]
         midSeq = chromSeq[midStart-1:midEnd]
         rightTSDSeq = chromSeq[tsd3Start-1:tsd3End]
         rightSeq = chromSeq[rightStart-1:rightEnd]
         
         res['emptySeq'] = leftSeq + leftTSDSeq + rightSeq
         emptyCoords = [tsd5Start,tsd5End]
         emptyCoords = [i - leftStart + 1 for i in emptyCoords]
         res['emptyTSDCoords'] = f'{emptyCoords[0]}-{emptyCoords[1]}'
         
         res['filledSeq'] = leftSeq + leftTSDSeq + midSeq + rightTSDSeq + rightSeq
         
         filledCoords = [tsd5Start,tsd5End,tsd3Start,tsd3End]
         filledCoords = [i - leftStart + 1 for i in filledCoords]
         tsdL = f'{filledCoords[0]}-{filledCoords[1]}'
         tsdR = f'{filledCoords[2]}-{filledCoords[3]}'
         res['filledTSDCoords'] = tsdL + ',' + tsdR
         
    else:
         leftEnd = tsd3Start - 1
         leftStart = leftEnd - flankSize + 1
         midStart = tsd3End + 1
         midEnd = tsd5Start - 1
         rightStart = tsd5End + 1
         rightEnd = rightStart + flankSize - 1
         
         leftSeq = chromSeq[leftStart-1:leftEnd]
         leftTSDSeq = chromSeq[tsd5Start-1:tsd5End]
         midSeq = chromSeq[midStart-1:midEnd]
         rightTSDSeq = chromSeq[tsd3Start-1:tsd3End]
         rightSeq = chromSeq[rightStart-1:rightEnd]
         
         res['emptySeq'] = leftSeq + rightTSDSeq + rightSeq
         res['emptySeq'] = refMEUtils.revcomp(res['emptySeq'])
         
         emptyCoords = [tsd3Start,tsd3End]
         emptyCoords = [i - leftStart + 1 for i in emptyCoords]
         emptyCoords = emptyCoords[::-1]
         # convert them to rc cords
         lenS = len(res['emptySeq'])
         emptyCoords = [ (lenS-i)+1  for i in emptyCoords]
         res['emptyTSDCoords'] = f'{emptyCoords[0]}-{emptyCoords[1]}'
         
         res['filledSeq'] = leftSeq + leftTSDSeq + midSeq + rightTSDSeq + rightSeq
         res['filledSeq'] = refMEUtils.revcomp(res['filledSeq'])
         
         filledCoords = [tsd3Start,tsd3End,tsd5Start,tsd5End]
         filledCoords = [i - leftStart + 1 for i in filledCoords]
                  
         filledCoords = filledCoords[::-1]
         # convert them to rc cords
         lenS = len(res['filledSeq'])
         filledCoords = [ (lenS-i)+1  for i in filledCoords]
         
         tsdL = f'{filledCoords[0]}-{filledCoords[1]}'
         tsdR = f'{filledCoords[2]}-{filledCoords[3]}'
         res['filledTSDCoords'] = tsdL + ',' + tsdR
         
    res['emptySeq'] = res['emptySeq'].upper()
    res['filledSeq'] = res['filledSeq'].upper()    
    return res 
#####################################################################        
# SETUP
parser = argparse.ArgumentParser(description='extract element sequecnes and empty sites')


parser.add_argument('--fasta', type=str,help='genome fata file',required=True)
parser.add_argument('--intable', type=str,help='input file of segments',required=True)
parser.add_argument('--outpre', type=str,help='prefix of output files',required=True)

args = parser.parse_args()
#####################################################################




print('input table:',args.intable)
print('output prefix:',args.outpre)


# read in sequences
genomeSeqs = refMEUtils.read_fasta_file_to_dict(args.fasta)
n = list(genomeSeqs.keys())
numSeqs = len(n)
print(f'read in {numSeqs} from fasta file')

outELEM_name =  args.outpre + '.elements.fa'
outELEM = open(outELEM_name,'w')
outALLELES_name = args.outpre + '.alleles.fa'
outALLELES = open(outALLELES_name,'w')

outTABLE_name = args.outpre + '.alleles.table'
outTABLE = open(outTABLE_name,'w')


numDid = 0

outTABLE.write('#locusID\tTSD_in_empty\tTSDs_in_filled\n')
inFile = open(args.intable,'r')
for line in inFile:
    ol = line
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
        continue

    numDid += 1
    if numDid % 1000 == 0:
        print(numDid)
        print(line)
        
    # need to process each element....    
    d = {}
    d['chrom'] = line[colToNames['chrom']]
    d['ori'] = line[colToNames['ori']]
    d['locusID'] = line[colToNames['locusID']]
    
    coords = line[colToNames['TSD-5prime-coords']]
    coords = coords.split('-')
    d['TSD-5prime-coords'] = [int(coords[0]),int(coords[1])]

    coords = line[colToNames['TSD-3prime-coords']]
    coords = coords.split('-')
    d['TSD-3prime-coords'] = [int(coords[0]),int(coords[1])]

    chromSeq = genomeSeqs[d['chrom']]['seq']
    
    
    # get element seq
    res = extract_element(d,chromSeq)    
    s = res['seq']
    s = refMEUtils.add_breaks_to_line(s,100)
    name = d['locusID']
    outELEM.write(f'>{name}\n{s}\n')
    
    # get filled site and empty site
    res = get_filled_empty(d,chromSeq,100)    
    
    #empty
    name = d['locusID'] + '_empty'
    s = res['emptySeq']
    s = refMEUtils.add_breaks_to_line(s,100)
    outALLELES.write(f'>{name}\n{s}\n')


    name = d['locusID'] + '_filled'
    s = res['filledSeq']
    s = refMEUtils.add_breaks_to_line(s,100)
    outALLELES.write(f'>{name}\n{s}\n')
    
    outTABLE.write('%s\t%s\t%s\n' % (d['locusID'],res['emptyTSDCoords'],res['filledTSDCoords']))
    




inFile.close()
outELEM.close()
outALLELES.close()
outTABLE.close()


