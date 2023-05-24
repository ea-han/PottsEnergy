import numpy as np
import io as io


class Energy:
    length = 0
    nVals = {}
    seqList = {}
    energy = 0

    def __init__(self, seqName, jName):
        #create energy object
        nFile = open(seqName, 'r')
        self.nVals = list(nFile.readline())
        self.length = len(self.nVals) - 1
        self.seqList = np.load(jName)
        self.energy = 0

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
        return 15
    
    def simplify(self):
        #simplifies nVals
        aList = ['R','H','K','D','E']
        bList = ['S','T','N','Q']
        cList = ['C','G','P']
        dList = ['A','V','I','L','M','F','Y','W']
        rawVals = list(self.nVals)
        for p in range(self.length):
            storChar = rawVals[p]
            if storChar in aList:
                self.nVals[p] = 'A'
            elif storChar in bList:
                self.nVals[p] = 'B'
            elif storChar in cList:
                self.nVals[p] = 'C'
            elif storChar in dList:
                self.nVals[p] = 'D'
            else:
                print("read incorrect char")
                return -1
    def getEnergy(self):
        ctr = 0
        for i in range(0,self.length):
            iAcid = self.nVals[i]
            for j in range (i+1,self.length):
                jAcid = self.nVals[j]
                couplingIdx = self.getCoupling(iAcid,jAcid)
                self.energy += self.seqList[ctr][couplingIdx]
                ctr+=1

    
    def genNpy(self, fileLen, fileName):
        newList = []
        ctr = 0
        for i in range(0,self.length):
            for j in range (i+1,self.length):
                if i < fileLen and j < fileLen:
                    newList.append(self.seqList[ctr])
                    ctr+=1
        np.save(fileName, newList, allow_pickle=False, fix_imports=False)
                



