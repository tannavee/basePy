def base_Capitalizer(seqList):
    uppercaseSeqList = []
    dnaChecker = True
    for nucleotide in seqList:
        if nucleotide.casefold() in ('a', 'c', 'g', 't'):
            uppercaseSeqList.append(nucleotide.upper())
        else:
            dnaChecker = False
    if dnaChecker == True:
        seqList = uppercaseSeqList
    else:
        return "Warning: This is not a DNA sequence."
    return seqList

def dna_compSeqMaker(dnaSeqList):
    three_to_fivecompSeqList = []
    for nucleotide in dnaSeqList:
        if nucleotide.casefold() == 't':
            three_to_fivecompSeqList.append('A')
                
        elif nucleotide.casefold() == 'a':
            three_to_fivecompSeqList.append('T')
                
        elif nucleotide.casefold() == 'c':
            three_to_fivecompSeqList.append('G')
           
        elif nucleotide.casefold() == 'g':
            three_to_fivecompSeqList.append('C')
        else:
            dnaChecker = False
            break           
    compSeqList = three_to_fivecompSeqList[::-1]
    return compSeqList
        
def dna_startCodonFinder(seqList, positionToStart):
    counter = positionToStart
    found = False
    while counter <= len(seqList)-1 and found is False:
        if counter + 2 <= len(seqList)-1 and seqList[counter : counter + 3] == ['A', 'T', 'G']:
            found = True
        else: 
            counter += 1      
    if found == True:
        return counter + 2
    else:
        return "No start codon was found"

def dna_orfFinder(dnaSeqList, positionToStart):
    thirdStartBase_position = dna_startCodonFinder(dnaSeqList, positionToStart)
    if type(thirdStartBase_position) == str:
        return thirdStartBase_position
    else: 
        counter = thirdStartBase_position + 1 
        found = False
        while counter <= len(dnaSeqList)-1 and found is False:
            if counter + 2 <= len(dnaSeqList)-1 and dnaSeqList[counter : counter + 3] in (['T', 'A', 'A'], ['T', 'A', 'G'], ['T', 'G', 'A']):
                found = True
            else:
                counter += 3
        if found == True:
            return [thirdStartBase_position - 2, counter + 2]
        else:
            return "No stop codon was found"

def dna_wholeSeqOrfFinder(dnaSeqList):
    positionToStart = 0
    orfsInSeq = {}
    orfLabel = 'ORF '
    orfCounter = 0
    while positionToStart <= len(dnaSeqList)-1:
        orf = dna_orfFinder(dnaSeqList, positionToStart)
        if type(orf) == list:
            orfCounter += 1
            orfsInSeq[orfLabel + str(orfCounter)] = orf
            positionToStart = orf[1] + 1
        else:
            if len(orfsInSeq) == 0:
                return orf
            else:
                break
    if len(orfsInSeq) > 1:
        orderedOrfs = {j: i for j, i in sorted(orfsInSeq.items(), key=lambda item: abs(item[1][0]-item[1][1]), reverse=True)}
        return orderedOrfs
    else:
        return orfsInSeq
    
    
def dna_orfPrinter(orfOutputDict):
    if type(orfOutputDict) != dict:
        return orfOutputDict 
    else:
        if len(orfOutputDict) == 1:
            printStatement = 'There is only one open reading frame in this sequence which spans from base ' + str(orfOutputDict['ORF 1'][0]) + ' to ' + str(orfOutputDict['ORF 1'][1]) + ' making it ' + str(orfOutputDict['ORF 1'][1] - orfOutputDict['ORF 1'][0]) + ' bases long'
            return printStatement
        else:
            printStatement = 'The following ' + str(len(orfOutputDict)) + ' open reading frames in the sequence are ordered from greatest to least in size: ' '\n'
            for key, (start, stop) in orfOutputDict.items():
                printStatement += '{}: {} bases long, spanning from bases {} to {}\n'.format(key, stop-start + 1, start, stop)
            return printStatement
    
def rna_codingToMRNA(dnaSeqList):
    mrnaList = []
    rnaChecker = True
    for nucleotide in dnaSeqList:
        if nucleotide.casefold() in ("a", "c", "g"):
            mrnaList.append(nucleotide.upper()) 
        elif nucleotide.casefold() == 't':
            mrnaList.append('U') 
        else:
            if nucleotide.casefold() not in  ("a", "t", "c", "g"):
                rnaChecker = False
                break              
    if rnaChecker == True: 
        return mrnaList
    else: 
        return "Warning: This is not a DNA sequence." 

def mRNA_rnaToAminoAcidSeq(list):
    counter = 0
    aminoAcidList = []
    while counter <= len(list) - 1:
        if counter == 0:
            aminoAcidList.extend(['Nter—', 'M'])
            counter += 3
        elif list[counter : counter + 3] in (['U', 'U', 'U'], ['U', 'U', 'C']):
            aminoAcidList.append('F')
            counter += 3
        elif list[counter : counter + 3] in (['U', 'U', 'A'], ['U', 'U', 'G'], ['C', 'U', 'U'], ['C', 'U', 'C'], ['C', 'U', 'A'], ['C', 'U', 'G']):
            aminoAcidList.append('L')
            counter += 3
        elif list[counter : counter + 3] in (['A', 'U', 'U'], ['A', 'U', 'C'], ['A', 'U', 'A']):
            aminoAcidList.append('I')
            counter += 3
        elif list[counter : counter + 3] == ['A', 'U', 'G']:
            aminoAcidList.append('M')
            counter += 3
        elif list[counter : counter + 3] in (['G', 'U', 'U'], ['G', 'U', 'C'], ['G', 'U', 'A'], ['G', 'U', 'G']):
            aminoAcidList.append('V')
            counter += 3
        elif list[counter : counter + 3] in (['U', 'C', 'U'], ['U', 'C', 'C'], ['U', 'C', 'A'], ['U', 'C', 'G'], ['A', 'G', 'U'], ['A', 'G', 'C']):
            aminoAcidList.append('S')
            counter += 3
        elif list[counter : counter + 3] in (['C', 'C', 'U'], ['C', 'C', 'C'], ['C', 'C', 'A'], ['C', 'C', 'G']):
            aminoAcidList.append('P')
            counter += 3
        elif list[counter : counter + 3] in (['A', 'C', 'U'], ['A', 'C', 'C'], ['A', 'C', 'A'], ['A', 'C', 'G']):
            aminoAcidList.append('T')
            counter += 3
        elif list[counter : counter + 3] in (['G', 'C', 'U'], ['G', 'C', 'C'], ['G', 'C', 'A'], ['G', 'C', 'G']):
            aminoAcidList.append('A')
            counter += 3
        elif list[counter : counter + 3] in (['U', 'A', 'U'], ['U', 'A', 'C']):
            aminoAcidList.append('Y')
            counter += 3
        elif list[counter : counter + 3] in (['C', 'A', 'U'], ['C', 'A', 'C']):
            aminoAcidList.append('H')
            counter += 3
        elif list[counter : counter + 3] in (['C', 'A', 'A'], ['C', 'A', 'G']):
            aminoAcidList.append('Q')
            counter += 3
        elif list[counter : counter + 3] in (['A', 'A', 'U'], ['A', 'A', 'C']):
            aminoAcidList.append('N')
            counter += 3
        elif list[counter : counter + 3] in (['A', 'A', 'A'], ['A', 'A', 'G']):
            aminoAcidList.append('K')
            counter += 3
        elif list[counter : counter + 3] in (['G', 'A', 'U'], ['G', 'A', 'C']):
            aminoAcidList.append('D')
            counter += 3
        elif list[counter : counter + 3] in (['G', 'A', 'A'], ['G', 'A', 'G']):
            aminoAcidList.append('E')
            counter += 3
        elif list[counter : counter + 3] in (['U', 'G', 'U'], ['U', 'G', 'C']):
            aminoAcidList.append('C')
            counter += 3
        elif list[counter : counter + 3] == ['U', 'G', 'G']:
            aminoAcidList.append('W')
            counter += 3
        elif list[counter : counter + 3] in (['C', 'G', 'U'], ['C', 'G', 'C'], ['C', 'G', 'A'], ['C', 'G', 'G'], ['A', 'G', 'A'], ['A', 'G', 'G']):
            aminoAcidList.append('R')
            counter += 3
        elif list[counter : counter + 3] in (['G', 'G', 'U'], ['G', 'G', 'C'], ['G', 'G', 'A'], ['G', 'G', 'G']):
            aminoAcidList.append('G')
            counter += 3            
        elif list[counter : counter + 3] in (['U', 'A', 'A'], ['U', 'A', 'G'], ['U', 'G', 'A']):
                aminoAcidList.append('—Cter')
                counter += 3
        else:
            aminoAcidList.append('ERROR')
            counter += 3
    return aminoAcidList     

def protein_printer(proteinOutputDict):
    if type(proteinOutputDict) != dict:
        return 'Warning: This is not the proper input.'
    else:
        if len(proteinOutputDict) == 1:
            printStatement = 'There is only one protein in the used sequence which is ' + str(len(proteinOutputDict[key]) - 2) + ' bases long'
            return printStatement
        else:
            printStatement = 'The following ' + str(len(proteinOutputDict)) + ' proteins in the sequence are ordered from greatest to least in size: ' '\n'
            for key, value in proteinOutputDict.items():
                printStatement += '{}: {} amino acids long\n'.format(key, len(value) - 2)
            return printStatement
