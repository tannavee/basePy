helix_chouFasmanDict = {
  'Ala': 1.29,
  'Cys': 1.11,
  'Leu': 1.3,
  'Met': 1.47,
  'Glu': 1.44,
  'Gln': 1.27,
  'His': 1.22,
  'Lys': 1.23,
  'Val': 0.91,
  'Ile': 0.97,
  'Phe': 1.07,
  'Tyr': 0.72,
  'Trp': 0.99,
  'Thr': 0.82,
  'Gly': 0.56,
  'Ser': 0.82,
  'Asp': 1.04,
  'Asn': 0.9,
  'Pro': 0.52,
  'Arg': 0.96,
  'A': 1.29,
  'C': 1.11,
  'L': 1.3,
  'M': 1.47,
  'E': 1.44,
  'Q': 1.27,
  'H': 1.22,
  'K': 1.23,
  'V': 0.91,
  'I': 0.97,
  'F': 1.07,
  'Y': 0.72,
  'W': 0.99,
  'T': 0.82,
  'G': 0.56,
  'S': 0.82,
  'D': 1.04,
  'N': 0.9,
  'P': 0.52,
  'R': 0.96
}

# converts all amino acids to their propensities
def pro_ProteintoHelixPropensity(aminoAcidSeq):
    propensityList = []
    for i in aminoAcidSeq:
        if i in ('Nter—', '—Cter'):
            label = i
        else: 
            aminoAcid = i
            propensityList.append(helix_chouFasmanDict[aminoAcid])
    return propensityList


# finds the start site for the alpha helix
def pro_initFinder(propensityList, positionToStart): #positionToStart):
    counter = positionToStart 
    found = False 
    while counter <= len(propensityList)-1 and found is False:
        # trying to find 4 amino acids in a row that are all > 1
        if counter + 3 <= len(propensityList)-1 and propensityList[counter] >= 1 and propensityList[counter + 1] >= 1 and propensityList[counter + 2] >= 1 and propensityList[counter + 3] >= 1:
            found = True 
        else:
            counter += 1
    if found == True:
        # will return a list with the start and end of the start site
        return [counter, counter + 3]
    else: 
        return "there is no initiation site"

# will expand the helix in the left direction, adds 2 to the inital 4, then begins adding 1, and taking the avg of the new + 3 at the end for > 1
def checkLeft(propensityList, positionToStart):
    startSite = pro_initFinder(propensityList, positionToStart)
    # if no initiation state found, return error
    if type(startSite) == str:
        return startSite 
    else:
        # add two AA
        expandedLeft = startSite
        if startSite[0] >= 2:
            expandedLeft[0] = expandedLeft[0] - 2
            while expandedLeft[0] >= 1:
                # adding 1 AA, and taking avg of it and 3 end AA for > 1
                if (propensityList[expandedLeft[0]-1] + propensityList[expandedLeft[0]] + propensityList[expandedLeft[0]+1] + propensityList[expandedLeft[0]+2])/4 >= 1:
                    expandedLeft[0] = expandedLeft[0] - 1
                else:
                    # will return the expanded list if avg != 1
                    return expandedLeft
            if expandedLeft[0] == 0:
                return expandedLeft
        else:
            return startSite

# expand to the right, and check that whole sequence avg > 1
def checkWhole(propensityList, positionToStart):
    # getting the expanded left section
    expandedLeft = checkLeft(propensityList, positionToStart)
    # if expanded left is error, return error
    if type(expandedLeft) == str:
        return expandedLeft
    else: 
        expandedRight = expandedLeft
        # adding 2 AA to the right side of the expanded left
        if (len(propensityList) - 1 - expandedRight[1]) >= 2:
            expandedRight[1] = expandedRight[1] + 2
            # adding 1 AA, and taking avg of it and 3 end AA for > 1
            while (len(propensityList) - 1 - expandedRight[1]) >= 1:
                if (propensityList[expandedRight[1]+1] + propensityList[expandedRight[1]] + propensityList[expandedRight[1]-1] + propensityList[expandedRight[1]-2])/4 >= 1:
                    expandedRight[1] = expandedRight[1] + 1
                else:
                    expandedWhole = expandedRight
                    # getting avg of whole seq, to make sure > 1
                    avg = sum(propensityList[expandedWhole[0]: expandedWhole[1] + 1]) / len(propensityList[expandedWhole[0]: expandedWhole[1] + 1])
                    if avg >= 1:
                        return expandedWhole
                    else:
                        # if not greater than one, return the position of the failed helix in tuple type
                        return tuple(expandedWhole)
            if (len(propensityList) - 1 - expandedRight[1]) == 0:
                expandedWhole = expandedRight
                avg = sum(propensityList[expandedWhole[0]: expandedWhole[1] + 1]) / len(propensityList[expandedWhole[0]: expandedWhole[1] + 1])
                if avg >= 1:
                    return expandedWhole
                else:
                    return tuple(expandedWhole)
        else:
            expandedWhole = expandedLeft
            # if size of list is 6, but start index is 1, and end start site is 4, this is too short 
            if expandedWhole[0] == 1 and expandedWhole[1] == 4 and len(propensityList) == 6:
                return "too short"
            else: 
                return expandedWhole

# function that returns a dictionary with postions of all the helices
def wholeSeqHelixFinder(propensityList):
    positionToStart = 0
    helicesInSeq = {}
    helixLabel = 'Helix '
    helixCounter = 0
    while positionToStart <= len(propensityList)-1:
        helix = checkWhole(propensityList, positionToStart)
        # checks to see if type list, if so will add to dictionary
        if type(helix) == list:
            helixCounter += 1
            helicesInSeq[helixLabel + str(helixCounter)] = helix
            positionToStart = helix[1] + 1
        # if tuple, will take the break position, and start looking for start site next position
        elif type(helix) == tuple:
            positionToStart = helix[1] + 1
        # is no helices in dectionary, return error msg
        else:
            if len(helicesInSeq) == 0:
                return 'no helices found'
            else:
                return helicesInSeq
    return helicesInSeq
    

        
