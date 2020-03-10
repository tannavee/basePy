import basePy as bp

def dna2protein(seqFile):
    # opening the txt file 
    fileContents = open(seqFile, 'r')
    # reading txt file, replaing newline characters with "", and blank spots with ""
    givenSeq = fileContents.read().replace('\n', '').replace(' ', '').replace('\t', '')
    # making the given sequence into a list
    givenSeqList = list(givenSeq)
    # making all bases uppercase, and checking that provided seq is valid
    givenSeqList = bp.base_Capitalizer(givenSeqList) 
    if type(givenSeqList) == str:
        return givenSeqList
    # if DNA sequence
    else:
        # using our library to create a complementary sequence list
        compSeqList = bp.dna_compSeqMaker(givenSeqList)
        # finding orfs for the given/comp seqs
        orfsInGivenSeq = bp.dna_wholeSeqOrfFinder(givenSeqList)
        orfsInCompSeq = bp.dna_wholeSeqOrfFinder(compSeqList)
        
        #if given doesn't have orf, and comp does
        if type(orfsInGivenSeq) == str and type(orfsInCompSeq) == dict:
            print("Results for the given sequence:", orfsInGivenSeq)
            print("Results for the complementary sequence:", bp.dna_orfPrinter(orfsInCompSeq))
            codingStrandList = compSeqList
            mRNApositions = {key.replace("ORF", "mRNA"): value for key, value in orfsInCompSeq.items()}

        
        # if given has orf and comp does
        elif type(orfsInCompSeq) == str and type(orfsInGivenSeq) == dict:
            print("Results for the complementary sequence:", orfsInCompSeq)
            print("Results for the given sequence:", bp.dna_orfPrinter(orfsInGivenSeq))
            codingStrandList = givenSeqList
            mRNApositions = {key.replace("ORF", "mRNA"): value for key, value in orfsInGivenSeq.items()}

            
        # if given and comp both have orfs
        elif type(orfsInCompSeq) == dict and type(orfsInGivenSeq) == dict:
            print("There are ORFS in the given and complementary sequences.")
            if (orfsInCompSeq['ORF 1'][1] - orfsInCompSeq['ORF 1'][0]) > (orfsInGivenSeq['ORF 1'][1] - orfsInGivenSeq['ORF 1'][0]):
                codingStrandList = compSeqList
                print('The longest ORF is in the complementary sequence, this program will use this sequence to get the mRNA.')
                print("Results for the complementary sequence:", bp.dna_orfPrinter(orfsInCompSeq))
                mRNApositions = mydictionary = {key.replace("ORF", "mRNA"): value for key, value in orfsInCompSeq.items()}

            else:
                codingStrandList = givenSeqList
                print('The longest ORF is in the given sequence, this program will use this sequence to get the mRNA.')
                print("Results for the given sequence:", bp.dna_orfPrinter(orfsInGivenSeq))
                mRNApositions = {key.replace("ORF", "mRNA"): value for key, value in orfsInGivenSeq.items()}

        # if neither given and comp have orfs
        else:
            return "Results for the given sequence: " + orfsInGivenSeq, "Results for the complementary sequence: " + orfsInCompSeq
    
    # replacing the 't' with 'u' in the coding strand
    rnaStrandList = bp.rna_codingToMRNA(codingStrandList)    
    
    # creating new dictionary where each mRNA key will have the whole sequence
    mRNAseqDict = {}
    for key, value in mRNApositions.items():
        mRNAseqDict[key] = rnaStrandList[mRNApositions[key][0]:mRNApositions[key][1] + 1]
    
    # converting the codons to amino acid
    proteinSeqDict = {}
    for key, value in mRNAseqDict.items():
        proteinSeqDict[key.replace("mRNA", "Protein")] = bp.mRNA_rnaToAminoAcidSeq(mRNAseqDict[key])
    print(bp.protein_printer(proteinSeqDict))
    print('This is the largest protein found: ' + ''.join(proteinSeqDict[list(proteinSeqDict.keys())[0]]))
    return proteinSeqDict
    
    
    
        
            

     
