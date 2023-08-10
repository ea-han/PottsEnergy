import os

class EnergyAverage:
    formattedMutation1 = ""
    formattedMutation2 = ""
    totalSeq = []
    consensus = []
    mutation1 = []
    mutation2 = []
    foundSeqs = []
    jFileArray = []
    reduxName = ""
    ##search list of sequences for mutations 1 and 2
    ##add found mutations to a list of sequences to calculate energy of
    ##reduce sequences and mutations with redux file
    ##calculate average of all sequences
    ##mutations: [index, source, [mut1, mut2,etc]]

    def __init__(self,totalSeq, consensus, mutation1, mutation2, jFileArray, reduxName):
        self.foundSeqs = []
        self.totalSeq = totalSeq
        self.reduxName = reduxName
        self.mutation1 = mutation1
        self.mutation2 = mutation2
        self.jFileArray = jFileArray
        self.consensus = self.reduceSequence(consensus)


        self.formattedMutation1 = self.arrToString(mutation1)
        self.formattedMutation2 = self.arrToString(mutation2)
        

    def arrToString(self,mutation):
        altMutations = '/'.join(mutation[2])
        newString = ''
        newString = newString + str(mutation[1])
        newString = newString + str(mutation[0]+1)
        newString = newString + altMutations
        return newString


    ## reduces mutations
    def reduceMutations(self): 
        ctr = 0
        with open(self.reduxName) as kFile:
            for line in kFile:
                if ctr == self.mutation1[0]:
                    line = line.strip()
                    insertArr = line.split()
                    aAlpha = list(insertArr[1])
                    bAlpha = list(insertArr[2])
                    cAlpha = list(insertArr[3])
                    dAlpha = list(insertArr[4])
                    
                    val1 = self.mutation1[1]
                    val2 = self.mutation1[2][0]
                    if val1 in aAlpha:
                        self.mutation1[1] = 'A'
                    elif val1 in bAlpha:
                        self.mutation1[1] = 'B'
                    elif val1 in cAlpha:
                        self.mutation1[1] = 'C'
                    elif val1 in dAlpha:
                        self.mutation1[1] = 'D'
                    else:
                        return -1
                    
                    if val2 in aAlpha:
                        self.mutation1[2] = 'A'
                    elif val2 in bAlpha:
                        self.mutation1[2] = 'B'
                    elif val2 in cAlpha:
                        self.mutation1[2] = 'C'
                    elif val2 in dAlpha:
                        self.mutation1[2] = 'D'
                    else:
                        return -1
                    
                    
                    
                elif ctr == self.mutation2[0]:
                    line = line.strip()
                    insertArr = line.split()
                    aAlpha = list(insertArr[1])
                    bAlpha = list(insertArr[2])
                    cAlpha = list(insertArr[3])
                    dAlpha = list(insertArr[4])
                    
                    val1 = self.mutation2[1]
                    val2 = self.mutation2[2][0]
                    if val1 in aAlpha:
                        self.mutation2[1] = 'A'
                    elif val1 in bAlpha:
                        self.mutation2[1] = 'B'
                    elif val1 in cAlpha:
                        self.mutation2[1] = 'C'
                    elif val1 in dAlpha:
                        self.mutation2[1] = 'D'
                    else:
                        return -1
                    
                    if val2 in aAlpha:
                        self.mutation2[2] = 'A'
                    elif val2 in bAlpha:
                        self.mutation2[2] = 'B'
                    elif val2 in cAlpha:
                        self.mutation2[2] = 'C'
                    elif val2 in dAlpha:
                        self.mutation2[2] = 'D'
                    else:
                        return -1
                ctr += 1
    
    ##reduces sequences
    def reduceSequence(self, seq):
        ctr = 0
        with open(self.reduxName) as kFile:
                returnSeq = []
                for line in kFile:
                    ##simlify
                    line = line.strip()
                    insertArr = line.split()
                    ctr = int(insertArr[0])
                    aAlpha = list(insertArr[1])
                    bAlpha = list(insertArr[2])
                    cAlpha = list(insertArr[3])
                    dAlpha = list(insertArr[4])
                    
                    val = seq[ctr]
                    if val in aAlpha:
                        returnSeq.append('A')
                    elif val in bAlpha:
                        returnSeq.append('B')
                    elif val in cAlpha:
                        returnSeq.append('C')
                    elif val in dAlpha:
                        returnSeq.append('D')
                    else:
                        print("`invalid char`:",val ,flush=True)
                        return -1
                    ctr+=1
        return returnSeq

    ##searches for sequences
    def searchSequences(self):
        ctr = 0
        for seq in self.totalSeq:
            if (seq[self.mutation1[0]] in self.mutation1[2]) and (seq[self.mutation2[0]] in self.mutation2[2]):
                self.foundSeqs.append(seq)
            ctr+=1

        if (len(self.foundSeqs) == 0):
            print("no sequences found...")
            return -1      

    def reduceFoundSeqs(self):
        for i in range(0, len(self.foundSeqs)):
            self.foundSeqs[i] = self.reduceSequence(self.foundSeqs[i])
            print(self.foundSeqs[i])

    def reduceTotalSeqs(self):
        for i in range(0, len(self.totalSeq)):
            self.totalSeq[i] = self.reduceSequence(self.totalSeq[i])
            print(self.totalSeq[i])

        
    def getCoupling(self, iVal,jVal):
        aVali = ord(iVal) - 65
        aValj = ord(jVal) - 65
        if (aVali < 0 or aVali > 3):
            print("invalid coupling")
            return
        if (aVali == 0 and aValj == 0):
            return 0
        if (aVali == 0 and aValj == 1):
            return 1
        if (aVali == 0 and aValj == 2):
            return 2
        if (aVali == 0 and aValj == 3):
            return 3
        if (aVali == 1 and aValj == 0):
            return 4
        if (aVali == 1 and aValj == 1):
            return 5
        if (aVali == 1 and aValj == 2):
            return 6
        if (aVali == 1 and aValj == 3):
            return 7
        if (aVali == 2 and aValj == 0):
            return 8
        if (aVali == 2 and aValj == 1):
            return 9
        if (aVali == 2 and aValj == 2):
            return 10
        if (aVali == 2 and aValj == 3):
            return 11
        if (aVali == 3 and aValj == 0):
            return 12
        if (aVali == 3 and aValj == 1):
            return 13
        if (aVali == 3 and aValj == 2):
            return 14
        if (aVali == 3 and aValj == 3):
            return 15
        return -1

    ##gets energy of a sequence
    def getEnergy(self, seq):
        ctr = 0
        energy = 0
        length = len(seq)

        for i in range(0,length):
            iAcid = seq[i]
            for j in range (i+1,length):
                jAcid = seq[j]
                couplingIdx = self.getCoupling(iAcid,jAcid)
                if couplingIdx < 0:
                    print("coupling failed")
                    return -1
                energy += self.jFileArray[ctr][couplingIdx]
                ctr+=1
        return energy

    ##undoes the mutation
    def mutateSequence1(self, seq,mutation):
        sequence = list(seq)
        sequence[mutation[0]] = mutation[1]
        return sequence

    ##undoes both mutations
    def mutateSequence2(self, seq):
        sequence = list(seq)
        sequence[self.mutation1[0]] = self.mutation1[1]
        sequence[self.mutation2[0]] = self.mutation2[1]
        return sequence
    


    def mutateSequence1Consensus(self, seq,mutation):
        sequence = list(seq)
        sequence[mutation[0]] = mutation[2]
        return sequence
    
    def mutateSequence2Consensus(self, seq):
        sequence = list(seq)
        sequence[self.mutation1[0]] = self.mutation1[2]
        sequence[self.mutation2[0]] = self.mutation2[2]
        return sequence

    ##Calculates energies and returns them in an array
    def calcEnergies(self, sequence): 
        defaultEnergy = self.getEnergy(self.mutateSequence2(sequence)) ##works
        m2Energy = self.getEnergy(self.mutateSequence1(sequence, self.mutation1))
        m1Energy = self.getEnergy(self.mutateSequence1(sequence, self.mutation2))
        m1m2Energy = self.getEnergy(sequence) ##works

        m1Delta = defaultEnergy - m1Energy 
        m2Delta = defaultEnergy - m2Energy
        m1m2Delta = defaultEnergy - m1m2Energy ##works
        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        m1Delta += deltaDeltaE
        m2Delta += deltaDeltaE
        return [deltaDeltaE,m1m2Delta,m1Delta,m2Delta]
    
    def calcEnergies2(self, sequence): 
        defaultEnergy = self.getEnergy(sequence) ##works
        m1Energy = self.getEnergy(self.mutateSequence1Consensus(sequence, self.mutation1))
        m2Energy = self.getEnergy(self.mutateSequence1Consensus(sequence, self.mutation2))
        m1m2Energy = self.getEnergy(self.mutateSequence2Consensus(sequence)) ##works

        m1Delta = defaultEnergy - m1Energy 
        m2Delta = defaultEnergy - m2Energy
        m1m2Delta = defaultEnergy - m1m2Energy 
        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        return [deltaDeltaE,m1m2Delta,m1Delta,m2Delta]

    def calcConsensus(self):
        defaultEnergy = self.getEnergy(self.consensus)
        m1Energy = self.getEnergy(self.mutateSequence1Consensus(self.consensus, self.mutation1))
        m2Energy = self.getEnergy(self.mutateSequence1Consensus(self.consensus, self.mutation2))
        m1m2Energy = self.getEnergy(self.mutateSequence2Consensus(self.consensus))

        m1Delta = defaultEnergy - m1Energy 
        m2Delta = defaultEnergy - m2Energy
        m1m2Delta = defaultEnergy - m1m2Energy 
        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)
        headerName = self.formattedMutation1 + '-' + self.formattedMutation2


        return [headerName, deltaDeltaE,m1m2Delta,m1Delta,m2Delta]


    def calcOne(self,seq,name):
        self.reduceMutations()
        seq = self.reduceSequence(seq)

        
        consensusArr = self.calcConsensus()
        ddConsensus = consensusArr[0]
        d1Consensus = consensusArr[1]
        d2Consensus = consensusArr[2]
        d12Consensus = consensusArr[3]   

        print("CONSENSUS: ")
        print("DDelta: ",ddConsensus)
        print("Delta12: ",d12Consensus)
        print("Delta1: ",d1Consensus)
        print("Delta2: ",d2Consensus, "\n")

        print("SEQUENCE: ", name, ", MUTATION: ",self.mutation1, self.mutation2)
        calcArray = self.calcEnergies2(seq)
        ddSum = calcArray[0]
        d1Sum = calcArray[1]
        d2Sum = calcArray[2]
        d12Sum = calcArray[3]
        print("DDelta: ",ddSum)
        print("Delta12: ",d12Sum)
        print("Delta1: ",d1Sum)
        print("Delta2: ",d2Sum)

        print('\n\n')


    def calcAverage(self):
        ddSum = 0
        d1Sum = 0
        d2Sum = 0
        d12Sum = 0
        if self.searchSequences() == -1:
            return
        self.reduceFoundSeqs()
        self.reduceMutations()
        len = 0
        
        for seq in self.foundSeqs:
            calcArray = self.calcEnergies(seq)
            ddSum += calcArray[0]
            d1Sum += calcArray[1]
            d2Sum += calcArray[2]
            d12Sum += calcArray[3]
            len+=1

        ddSum = ddSum/len
        d1Sum = d1Sum/len
        d2Sum = d2Sum/len
        d12Sum = d12Sum/len

        consensusArr = self.calcConsensus()
        ddConsensus = consensusArr[0]
        d1Consensus = consensusArr[1]
        d2Consensus = consensusArr[2]
        d12Consensus = consensusArr[3]     

        print("AVERAGE: ")
        print("number of ",self.mutation1, self.mutation2, ": ", len)
        print("DDelta: ",ddSum)
        print("Delta12: ",d12Sum)
        print("Delta1: ",d1Sum)
        print("Delta2: ",d2Sum)
        print("\nCONSENSUS: ")
        print("DDelta: ",ddConsensus)
        print("Delta12: ",d12Consensus)
        print("Delta1: ",d1Consensus)
        print("Delta2: ",d2Consensus, "\n\n")

    def calcOneToOutput(self,seq,name):
        self.reduceMutations()
        seq = self.reduceSequence(seq)

        calcArray = self.calcEnergies2(seq)
        ddSum = calcArray[0]
        d12Sum = calcArray[1]
        d1Sum = calcArray[2]
        d2Sum = calcArray[3]
        headerName = self.formattedMutation1 + '-' + self.formattedMutation2


        returnArr = [headerName,name,ddSum,d12Sum,d1Sum,d2Sum]

        return returnArr





        

        