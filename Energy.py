import numpy as np
import io as io
import sys

class Energy:
    nVals = {} ##should not be modified, unsimplified sequence.
    seqList = {} ##numpy file
    posTuples = [] ##list of mutations
    outputFilename = ""
    reduxName = ""


    def __init__(self, seqName, jName, mutList, reduxName, outputFilename):
        #create energy object
        nFile = open(seqName, 'r')
        ##loads nVals with input sequence
        self.nVals = list(nFile.readline())
        self.nVals.pop(len(self.nVals)-1)
        self.seqList = np.load(jName)
        self.posTuples = []

        ##loads mutation list array
        with open(mutList) as kFile:
            for line in kFile:
                line = line.strip()
                insertArr = line.split(',')
                insertArr[1] = int(insertArr[1])
                self.posTuples.append(insertArr)
        self.outputFilename = outputFilename
        self.reduxName = reduxName
        


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

    ##takes a list of mutations, gets modified sequence and gets the energy of that particular sequence+set of mutations
    def getEnergy(self, seq, mutationList):
        ctr = 0
        energy = 0
        inputList = self.initMutated(seq, mutationList)
        length = len(inputList)

        for i in range(0,length):
            iAcid = inputList[i]
            for j in range (i+1,length):
                jAcid = inputList[j]
                couplingIdx = self.getCoupling(iAcid,jAcid)
                if couplingIdx < 0:
                    print("coupling failed")
                    return -1
                energy += self.seqList[ctr][couplingIdx]
                ctr+=1
        return energy
    
    ##takes a sequence and reducer as input, returns reduced sequence. NOTE: aAlpha still has leading '-'...
    def reduce(self, inSeq, reduxName):
        returnSeq = list(inSeq)
        with open(reduxName) as kFile:
            for line in kFile:
                aAlpha = []
                bAlpha = []
                cAlpha = []
                dAlpha = []

                line = line.strip()
                insertArr = line.split()
                ctr = int(insertArr[0])
                aAlpha = list(insertArr[1])
                bAlpha = list(insertArr[2])
                cAlpha = list(insertArr[3])
                dAlpha = list(insertArr[4])
                val = inSeq[ctr]

                if val in aAlpha:
                    returnSeq[ctr] = 'A'
                elif val in bAlpha:
                    returnSeq[ctr] = 'B'
                elif val in cAlpha:
                    returnSeq[ctr] = 'C'
                elif val in dAlpha:
                    returnSeq[ctr] = 'D'
                else:
                    return -1
                ctr+=1
        return returnSeq


    ## takes list of mutations and returns a mutated sequence list
    def initMutated(self, inSeq, mutList):
        newList = list(inSeq)
        for pos in mutList:
            newList[(pos[1])-1] = pos[2]

        newList = self.reduce(newList, self.reduxName)
        return newList    
    
    ##takes average deltadeltaE of a particular pair of mutations.
    ##seqList[][] - array of sequence arrays
    def mutateAverageOf1(self,seqList):
        ##for every sequence, 
        ##if seq[m1] = k1 and seq[m2] = k2
        sum = 0
        mutList = self.posTuples
        m1Idx = mutList[0][1]
        m2Idx = mutList[1][1]
        wildTypeEnergy = self.getEnergy(self.nVals, [])
        
        ctr = 0
        for seq in seqList:
            if seq[m1Idx-1] == mutList[0][2] and seq[m2Idx-1] == mutList[1][2]: ##if the sequence has both mutations
                ##get energy of seq
                ##get energy of wildtype
                ##calculate delta e of seq and wildtype
                seqEnergy = self.getEnergy(seq, [])
                deltaE = wildTypeEnergy - seqEnergy
                sum += deltaE
                ctr += 1

        if (ctr == 0):
            print("No sequences found.")
            return
        
        print("Average for ", mutList, " is: ", sum/ctr)

    ##Flips mutations in mutList
    def reverseMutList(self, mutList):
        newList = list(mutList)
        for pos in newList:
            temp = pos[0]
            pos[0] = pos[2]
            pos[2] = temp
        return newList


    def mutateAverageOf2(self,seqList):
        ## (wildtype - energy(sequence)) - (wildtype - energy(sequence without m2)) - (wildtype-energy(sequence without m1))
        sum = 0
        ctr = 0
        mutList = self.posTuples
        m1Idx = mutList[0][1]
        m2Idx = mutList[1][1]
        print(mutList)
        print(m1Idx)
        print(m2Idx)
        for seq in seqList:
            print(seq)
            if seq[m1Idx - 1] == mutList[0][2] and seq[m2Idx - 1] == mutList[1][2]: ##if the sequence has both mutations
                sum += self.mutateReturnReversed(seq, mutList)
                ctr += 1

        if (ctr == 0):
            print("No sequences found.")
            return
        return sum/ctr
    
    def mutateAverageOf3(self,seqList):
        ##Search for sequences that can possible mutate to m1m2, calculate deltadeltaE
        sum = 0
        ctr = 0
        mutList = self.posTuples
        m1Idx = mutList[0][1]
        m2Idx = mutList[1][1]
        for seq in seqList:
            if seq[m1Idx - 1] == mutList[0][0] and seq[m2Idx - 1] == mutList[1][0]: ##if the sequence can have both mutations
                sum += self.mutateReturn(seq, mutList)
                ctr += 1

        if (ctr == 0):
            print("No sequences found.")
            return
        return sum/ctr

            

    def mutateReturnReversed(self, sequence, mutationList): ##1 based indexing, pos and mutation
        reversedList = self.reverseMutList(mutationList)
        revertm1MutList = [None] *1
        revertm2MutList = [None] *1
        
        revertm1MutList[0] = reversedList[1]
        revertm2MutList[0] = reversedList[0]
        revertNoneList = []

        wildTypeEnergy = self.getEnergy(self.nVals, [])
        m1Energy = self.getEnergy(sequence, revertm2MutList)
        m2Energy = self.getEnergy(sequence, revertm1MutList)
        m1m2Energy = self.getEnergy(sequence, revertNoneList)

        m1Delta = wildTypeEnergy - m1Energy
        m2Delta = wildTypeEnergy - m2Energy
        m1m2Delta = wildTypeEnergy - m1m2Energy

        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        return deltaDeltaE

    ##calculates deltadeltae and returns 
    def mutateReturn(self, sequence, mutationList): ##1 based indexing, pos and mutation
        m1MutList = [None] *1
        m2MutList = [None] *1
        m1m2MutList = mutationList
        m1MutList[0] = mutationList[0]
        m2MutList[0] = mutationList[1]

        defaultEnergy = self.getEnergy(sequence, [])
        m1Energy = self.getEnergy(sequence, m1MutList)
        m2Energy = self.getEnergy(sequence, m2MutList)
        m1m2Energy = self.getEnergy(sequence, m1m2MutList)

        m1Delta = defaultEnergy - m1Energy
        m2Delta = defaultEnergy - m2Energy
        m1m2Delta = defaultEnergy - m1m2Energy

        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        return deltaDeltaE


    def mutate(self): ##1 based indexing, pos and mutation
        m1MutList = [None] *1
        m2MutList = [None] *1
        m1m2MutList = self.posTuples
        m1MutList[0] = self.posTuples[0]
        m2MutList[0] = self.posTuples[1]

        wildTypeEnergy = self.getEnergy(self.nVals, [])
        m1Energy = self.getEnergy(self.nVals, m1MutList)
        m2Energy = self.getEnergy(self.nVals, m2MutList)
        m1m2Energy = self.getEnergy(self.nVals, m1m2MutList)

        m1Delta = wildTypeEnergy - m1Energy
        m2Delta = wildTypeEnergy - m2Energy
        m1m2Delta = wildTypeEnergy - m1m2Energy

        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        stdoutStore = sys.stdout

        sys.stdout = open(self.outputFilename, 'w')
        print("Mutations: ", self.posTuples)
        print("Wildtype Energy: ", wildTypeEnergy)
        print("m1 Energy: ", m1Energy)
        print("m2 Energy: ", m2Energy)
        print("m1m2 Energy: ", m1m2Energy)
        print("DeltaE m1m2: ", m1m2Delta)
        print("DeltaE m1: ", m1Delta)
        print("DeltaE m2: ", m2Delta)
        print("\nDeltaDelta E: ", deltaDeltaE)
        sys.stdout.flush()

        sys.stdout = stdoutStore




