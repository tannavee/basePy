def dna_compSeqMaker(dnaSeqList):
    three_to_fivecompSeqList = []
    dnaChecker = True
    for nucleotide in dnaSeqList:
        if nucleotide.casefold() == 't':
            if nucleotide.isupper():
                three_to_fivecompSeqList.append('A')
            else:
                three_to_fivecompSeqList.append('a')
                
        elif nucleotide.casefold() == 'a':
            if nucleotide.isupper():
                three_to_fivecompSeqList.append('T')
            else:
                three_to_fivecompSeqList.append('t')
                
        elif nucleotide.casefold() == 'c':
            if nucleotide.isupper():
                three_to_fivecompSeqList.append('G')
            else:
                three_to_fivecompSeqList.append('g')
           
        elif nucleotide.casefold() == 'g':
            if nucleotide.isupper():
                three_to_fivecompSeqList.append('C')
            else:
                three_to_fivecompSeqList.append('c')
        else:
            dnaChecker = False
            break           
    compSeqList = three_to_fivecompSeqList[::-1]
    if dnaChecker == True: 
        return compSeqList
    else: 
        return "Warning: This is not a DNA sequence."

def dna_startCodonFinder(seqList, positionToStart):
    counter = positionToStart
    firstStart = False
    secondStart = False
    thirdStart = False
    while counter <= len(seqList)-1 and firstStart is False and secondStart is False and thirdStart is False:
        if counter + 2 <= len(seqList)-1 and seqList[counter].casefold() == 'a' and seqList[counter + 1].casefold() == 't' and seqList[counter + 2].casefold() == 'g':
            firstStart, secondStart, thirdStart = True, True, True 
        else: 
            counter += 1      
    if thirdStart == True:
        return counter + 2
    else:
        return "No start codon was found"

def dna_orfFinder(dnaSeqList, positionToStart):
    thirdStartBase_position = dna_startCodonFinder(dnaSeqList, positionToStart)
    if thirdStartBase_position == "No start codon was found":
        return thirdStartBase_position
    else: 
        counter = thirdStartBase_position + 1 
        firstStop = False
        secondStop = False
        thirdStop = False
        while counter <= len(dnaSeqList)-1 and firstStop is False and secondStop is False and thirdStop is False:
            if counter + 2 <= len(dnaSeqList)-1 and dnaSeqList[counter].casefold() == 't' and dnaSeqList[counter + 1].casefold() == 'a' and dnaSeqList[counter + 2].casefold() == 'a':
                firstStop, secondStop, thirdStop = True, True, True
            elif counter + 2 <= len(dnaSeqList)-1 and dnaSeqList[counter].casefold() == 't' and dnaSeqList[counter + 1].casefold() == 'a' and dnaSeqList[counter + 2].casefold() == 'g':
                firstStop, secondStop, thirdStop = True, True, True
            elif counter + 2 <= len(dnaSeqList)-1 and dnaSeqList[counter].casefold() == 't' and dnaSeqList[counter + 1].casefold() == 'g' and dnaSeqList[counter + 2].casefold() == 'a':
                firstStop, secondStop, thirdStop = True, True, True
            else:
                counter += 3
        if thirdStop == True:
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
    
    
def dna_orfPrinter(orfOutput):
    if type(orfOutput) != dict:
        return orfOutput 
    else:
        if len(orfOutput) == 1:
            printStatement = 'There is only one open reading frame in this sequence which spans from base ' + str(orfOutput['ORF 1'][0]) + ' to ' + str(orfOutput['ORF 1'][1]) + ' making it ' + str(orfOutput['ORF 1'][1] - orfOutput['ORF 1'][0]) + ' bases long'
            return printStatement
        else:
            printStatement = 'The following ' + str(len(orfOutput)) + ' open reading frames in the sequence are ordered from greatest to least in size: ' '\n'
            for key, (start, stop) in orfOutput.items():
                printStatement += '{}: {} bases long, spanning from bases {} to {}\n'.format(key, stop-start, start, stop)
            return printStatement
    

def rna_codingToMRNA(dnaSeqList):
    mrnaList = []
    rnaChecker = True
    for nucleotide in dnaSeqList:
        if nucleotide.casefold() in  ("a", "c", "g"):
            mrnaList.append(nucleotide) 
        elif nucleotide.casefold() == 't':
            if nucleotide.isupper():
                mrnaList.append('U')
            else:
                mrnaList.append('u') 
        else:
            if nucleotide.casefold() not in  ("a", "t", "c", "g"):
                rnaChecker = False
                break              
    if rnaChecker == True: 
        return mrnaList
    else: 
        return "Warning: This is not a DNA sequence."    


    





# tests cases

noStart = ['c', 'a', 't', 'c', 'a', 't', 'c', 'a', 't', 'c', 'a', 't', 'c', 'a', 't', 'c', 'a', 't', 'c', 'a', 't']
noStop = ['c', 'a', 'a', 't', 'g', 'c', 'a', 'a',  'c', 'a', 'a',  'c', 'a', 'a',  'c', 'a', 'a',  'c', 'a', 'a',  'c', 'a', 'a']
manyOrfsNoExtra = ['c', 'a', 't', 'g', 'c', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'a', 't', 'g', 'c', 'a', 'a', 'c', 'a', 'a', 't', 'g', 'a']
orfsWithExtra = ['c', 'a', 't', 'g', 'c', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'a', 't', 'g', 'c', 'a', 'a', 'c', 'a', 'a', 't', 'g', 'a', 'g', 'g', 'g', 'g', 'g']
startNoStart = ['c', 'a', 't', 'g', 'c', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'a', 't', 'g', 'c', 'a', 'a', 'c', 'a', 'a', 't', 'g', 'a', 'a', 't', 'g', 'a', 't']
startAtEnd = ['c', 'a', 't', 'g', 'c', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'a', 't', 'g', 'c', 'a', 'a', 'c', 'a', 'a', 't', 'g', 'a', 'a', 't', 'g']
startNextStop = ['a', 't', 'g', 't', 'a', 'a']
oneOrfwithStarts = ['a', 't', 'g', 'a', 't', 'g', 'a', 't', 'g', 't', 'a', 'a']
manyStartNextStop = ['a', 't', 'g', 't', 'a', 'a', 'a', 't', 'g', 't', 'a', 'a']
manyStartNextStopwithFiller = ['a', 't', 'a', 't', 'g', 't', 'a', 'a', 'a', 't', 'a', 't', 'a', 't', 'g', 't', 'a', 'a']

        
        