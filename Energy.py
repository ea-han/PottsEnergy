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

        self.posTuples = self.reduceMutationList(self.posTuples, reduxName)
        self.nVals = self.reduce(self.nVals)


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
    def reduceMutationList(self, mutList, reduxName):
        returnList = list(mutList)
        idx1 = mutList[0][1]-1
        idx2 = mutList[1][1]-1

        with open(reduxName) as kFile:

                # lines to read
            line_numbers = [idx1, idx2]
            # To store lines
            lines = []
            for i, line in enumerate(kFile):
                # read lines
                if i in line_numbers:
                    lines.append(line.strip())
                elif i > idx2:
                    # don't read after line 7 to save time
                    break
            ctr = 0
            for line in lines:
                insertArr = line.split()
                aAlpha = list(insertArr[1])
                bAlpha = list(insertArr[2])
                cAlpha = list(insertArr[3])
                dAlpha = list(insertArr[4])
                preMut = mutList[ctr][0]
                postMut = mutList[ctr][2]

                if preMut in aAlpha:
                    returnList[ctr][0] = 'A'
                elif preMut in bAlpha:
                    returnList[ctr][0] = 'B'
                elif preMut in cAlpha:
                    returnList[ctr][0] = 'C'
                elif preMut in dAlpha:
                    returnList[ctr][0] = 'D'
                else:
                    return -1
                
                if postMut in aAlpha:
                    returnList[ctr][2] = 'A'
                elif postMut in bAlpha:
                    returnList[ctr][2] = 'B'
                elif postMut in cAlpha:
                    returnList[ctr][2] = 'C'
                elif postMut in dAlpha:
                    returnList[ctr][2] = 'D'
                else:
                    return -1
                ctr+=1
         
        return returnList
    

    def reduce(self, inSeq): ##takes a sequence and reducer as input, returns reduced sequence. NOTE: aAlpha still has leading '-'...
        returnSeq = list(inSeq)
        with open(self.reduxName) as kFile:
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
        return newList     
    
    ## reverses mutation lists
    def reverse(self,inputList):
        mutList = list(inputList)
        for pos in mutList:
            temp = pos[0]
            pos[0] = pos[2]
            pos[2] = temp
        return mutList
    

    def mutateAverageOf3(self,seqList):
        ##Search for sequences that have mutated
        delta2sum = 0
        m1Sum = 0
        m2Sum = 0
        m1m2Sum = 0
        ctr = 0
        mutList = self.posTuples
        print(mutList)
        m1Idx = mutList[0][1]
        m2Idx = mutList[1][1]

        for seq in seqList:
            seq = self.reduce(seq)
            if seq[m1Idx - 1] == mutList[0][2] and seq[m2Idx - 1] == mutList[1][2]: ##if the sequence has both mutations
                arr = self.mutateReturn(seq,self.reverse(mutList))
                delta2sum += arr[0]
                m1Sum += arr[1]
                m2Sum += arr[2]
                m1m2Sum += arr[3]
                ctr += 1
                
        print(ctr)
        if (ctr == 0):
            print("No sequences found.")
            return

        sys.stdout = open(self.outputFilename, 'w')
        print("Mutations: ", mutList)
        print("Delta-Delta E: ", delta2sum/ctr)
        print("Delta-m1m2: ", m1m2Sum/ctr)
        print("Delta-m1: ", m1Sum/ctr)
        print("Delta-m2: ", m2Sum/ctr)
        sys.stdout.flush()


    ##calculates deltadeltae and returns 
    def mutateReturn(self, sequence, mutationList): ##1 based indexing, pos and mutation
        m1MutList = [None] *1
        m2MutList = [None] *1
        m1m2MutList = mutationList
        m1MutList[0] = mutationList[0]
        m2MutList[0] = mutationList[1]

        defaultEnergy = self.getEnergy(sequence, m1m2MutList)
        m1Energy = self.getEnergy(sequence, m2MutList)
        m2Energy = self.getEnergy(sequence, m1MutList)
        m1m2Energy = self.getEnergy(sequence, [])

        m1Delta = defaultEnergy - m1Energy
        m2Delta = defaultEnergy - m2Energy
        m1m2Delta = defaultEnergy - m1m2Energy

        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        return [deltaDeltaE,m1Delta,m2Delta,m1m2Delta]


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




