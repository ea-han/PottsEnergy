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
                insertArr[0] = int(insertArr[0])
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
    def getEnergy(self, mutationList):
        ctr = 0
        energy = 0
        inputList = self.initMutations(mutationList)
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


    ## takes list of mutations and returns a mutated sequence list NOTE: commented out until I get files from avik
    def initMutations(self, mutList):
        newList = list(self.nVals)
        for pos in mutList:
            newList[(pos[0] - 1)] = pos[1]

        ##newList = self.reduce(newList, self.reduxName)
        return newList    

    def mutate(self): ##1 based indexing, pos and mutation
        m1MutList = [None] *1
        m2MutList = [None] *1
        m1m2MutList = self.posTuples
        m1MutList[0] = self.posTuples[0]
        m2MutList[0] = self.posTuples[1]

        wildTypeEnergy = self.getEnergy([])
        m1Energy = self.getEnergy(m1MutList)
        m2Energy = self.getEnergy(m2MutList)
        m1m2Energy = self.getEnergy(m1m2MutList)

        m1Delta = wildTypeEnergy - m1Energy
        m2Delta = wildTypeEnergy - m2Energy
        m1m2Delta = wildTypeEnergy - m1m2Energy

        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        stdoutStore = sys.stdout

        sys.stdout = open(self.outputFilename, 'w')
        print("Mutations: ", self.posTuples)
        print("Wildtype Energy: ", wildTypeEnergy)
        print("m1 Energy: ", wildTypeEnergy)
        print("m2 Energy: ", wildTypeEnergy)
        print("m1m2 Energy: ", wildTypeEnergy)
        print("DeltaE m1m2: ", m1m2Delta)
        print("DeltaE m1: ", m1Delta)
        print("DeltaE m2: ", m2Delta)
        print("\nDeltaDelta E: ", deltaDeltaE)
        sys.stdout.flush()

        sys.stdout = stdoutStore




