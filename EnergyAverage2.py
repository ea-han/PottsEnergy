import os
import pdb
import json
import pandas as pd
import csv

## TODO: 
# -get mutations/sequence pair where they are flipped to consensus
# -get unreduced mutations for csv file output
# -(optional) write two comparison functions, one for DistA and DistB
# -test code

class EnergyAverage2:
    seqStack = []
    consensus = []
    jFileArray = []
    reduxName = ""
    mutationList = []

    ##[valStack, valStack, valStack,...]
    importedStacks = []
    
    INTERACTION1 = "Rescue"
    INTERACTION2 = "Compensatory"
    INTERACTION3 = "Antagonistic"

    ##search list of sequences for mutations 1 and 2
    ##add found mutations to a list of sequences to calculate energy of
    ##reduce sequences and mutations with redux file
    ##calculate average of all sequences
    ##mutations: [index, source, [mut1, mut2,etc]],[index2, source2, [mut1, mut2,etc]],

    def __init__(self,seqStack, consensus, mutationList, jFileArray, reduxName):
        self.foundSeqs = []
        self.seqStack = seqStack
        self.reduxName = reduxName
        self.jFileArray = jFileArray
        self.consensus = self.reduceSequence(consensus)
        self.mutationList = mutationList

    def arrToString(self,mutation):
        altMutations = '/'.join(mutation[2])
        newString = ''
        newString = newString + str(mutation[1])
        newString = newString + str(mutation[0]+1)
        newString = newString + altMutations
        return newString
    
    def mutPairToString(self,mutationPair):
        string1 = self.arrToString(mutationPair[0])
        string2 = self.arrToString(mutationPair[1])
        return string1 + "-" + string2

    ##  reduces sequence stack to json file
    def preReduceToJson(self):
        newFullStack = []
        for seq in self.seqStack:
            newFullStack.append(self.reduceSequence(seq))
            with open("out/reduced.in.fullseq.json", "w") as outfile:
                json.dump(newFullStack, outfile)

    ##loads stack from json file
    def loadJson(self):
        with open("out/reduced.in.fullseq.json") as inputFile:
            self.seqStack = json.load(inputFile)

    def reduceAllMutations(self):
        for mutPair in self.mutationList:
            mutPair = self.reduceMutations(mutPair)
        #breakpoint()

    ## reduces mutations
    def reduceMutations(self,mutationPair): 
        ctr = 0
        mutation1 = mutationPair[0]
        mutation2 = mutationPair[1]
        with open(self.reduxName) as kFile:
            for line in kFile:
                if ctr == mutation1[0]:
                    line = line.strip()
                    insertArr = line.split()
                    aAlpha = list(insertArr[1])
                    bAlpha = list(insertArr[2])
                    cAlpha = list(insertArr[3])
                    dAlpha = list(insertArr[4])
                    
                    val1 = mutation1[1]
                    val2 = mutation1[2][0]
                    if val1 in aAlpha:
                        mutation1[1] = 'A'
                    elif val1 in bAlpha:
                        mutation1[1] = 'B'
                    elif val1 in cAlpha:
                        mutation1[1] = 'C'
                    elif val1 in dAlpha:
                        mutation1[1] = 'D'
                    else:
                        print("invalid mutation: " + mutation1)
                        return -1
                    
                    if val2 in aAlpha:
                        mutation1[2] = 'A'
                    elif val2 in bAlpha:
                        mutation1[2] = 'B'
                    elif val2 in cAlpha:
                        mutation1[2] = 'C'
                    elif val2 in dAlpha:
                        mutation1[2] = 'D'
                    else:
                        print("invalid mutation: " + mutation1)
                        return -1
                    
                    
                    
                elif ctr == mutation2[0]:
                    line = line.strip()
                    insertArr = line.split()
                    aAlpha = list(insertArr[1])
                    bAlpha = list(insertArr[2])
                    cAlpha = list(insertArr[3])
                    dAlpha = list(insertArr[4])
                    
                    val1 = mutation2[1]
                    val2 = mutation2[2][0]
                    if val1 in aAlpha:
                        mutation2[1] = 'A'
                    elif val1 in bAlpha:
                        mutation2[1] = 'B'
                    elif val1 in cAlpha:
                        mutation2[1] = 'C'
                    elif val1 in dAlpha:
                        mutation2[1] = 'D'
                    else:
                        print("invalid mutation: " + mutation2)
                        return -1
                    
                    if val2 in aAlpha:
                        mutation2[2] = 'A'
                    elif val2 in bAlpha:
                        mutation2[2] = 'B'
                    elif val2 in cAlpha:
                        mutation2[2] = 'C'
                    elif val2 in dAlpha:
                        mutation2[2] = 'D'
                    else:
                        print("invalid mutation: " + mutation2)
                        return -1
                ctr += 1

            
        
        return [mutation1,mutation2]
    
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

    def reduceStack(self):
        for i in range(0, len(self.seqStack)):
            self.seqStack[i] = self.reduceSequence(self.seqStack[i])

        
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

    def mutateSequence1Consensus(self, seq,mutation):
        sequence = list(seq)
        sequence[mutation[0]] = mutation[2]
        return sequence
    
    def mutateSequence2Consensus(self, seq,mutation1,mutation2):
        sequence = list(seq)
        sequence[mutation1[0]] = mutation1[2]
        sequence[mutation2[0]] = mutation2[2]
        return sequence
    
    ##takes reduced sequence andmutation and calcultes energies. Need to screen for valid mutation/sequences.
    #runs? need to check avlues
    def calcEnergySequence(self,seq,name,mutationArray):
        
        defaultEnergy = self.getEnergy(seq) ##works
        m1Energy = self.getEnergy(self.mutateSequence1Consensus(seq, mutationArray[0]))
        m2Energy = self.getEnergy(self.mutateSequence1Consensus(seq, mutationArray[1]))
        m1m2Energy = self.getEnergy(self.mutateSequence2Consensus(seq,mutationArray[0],mutationArray[1])) ##works

        m1Delta = defaultEnergy - m1Energy 
        m2Delta = defaultEnergy - m2Energy
        m1m2Delta = defaultEnergy - m1m2Energy 
        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        return [mutationArray,name,deltaDeltaE,m1m2Delta,m1Delta,m2Delta]
    
    ##computes all the sequence values for a given mutation
    ##TODO rewrite in C and call cdef from cython
    ##      -checkIfEligible, calcEnergySequence
    def computeSeqsForMutation(self,mutationPair):
        outputArr = []
        ##mutationpair, seq index, delta delta, d12, d1, d2
        index = 0
        for seq in self.seqStack:
            if self.checkIfEligible(seq, mutationPair) == True:
                outputArr.append(self.calcEnergySequence(seq,index,mutationPair))
            index+=1
        ##breakpoint() ##computed 1 seqstack
        return outputArr


    ##for each mutation in mutationlist, compute all values. return list of lists of values
    def computeValueStacks(self):
        ##reduces stack NOTE: now using imported json file
        ##self.reduceStack()
        ##output Format: [[mutation,computed stack values],[mutation,computed stack values]]
        computedArr = []
        outputArr = []

        #this runs
        for mutPair in self.mutationList:
            computedArr = self.computeSeqsForMutation(mutPair)
            ##add stack to outputarr
            outputArr.append([mutPair,computedArr])

            ##write to json
            ##os.mkdir("out/"+self.mutPairToString(mutPair))

            with open("out/"+self.mutPairToString(mutPair)+"/valueStack.json", "w") as outfile1:
                json.dump(outputArr, outfile1)

        print("computed Value Stacks")
        return outputArr
    

    
    def getHDNum(self,seq):
        seqMut = ''
        consMut = ''
        mutList = []
        for i in range(0,len(self.consensus)):
            seqMut = seq[i] 
            consMut = self.consensus[i]
            if seqMut != consMut:
                addString = consMut+str(i)+seqMut
                mutList.append(addString)

        return mutList

    ##works on consensus?
    def checkIfEligible(self, seq, mutPair):
        mut1 = mutPair[0]
        mut2 = mutPair[1]
        index1 = mut1[0]
        index2 = mut2[0]
        preMut1 = mut1[1]
        preMut2 = mut2[1]

        if seq[index1] == preMut1 and seq[index2] == preMut2:
            return True
        else:
            return False



    ##takes in SINGLE MUTATION valuestack,compares each sequence DDE with consensus and classifies.
    ##Sorts and returns all values in valuestack
    ##returns []
    ##REWRITE
    def getInteractionsAndDifference(self,mutationData):
        mutPair = mutationData[0]
        valueStack = mutationData[1]

        rescue = []
        compensate = []
        antag = []

        ##pdb.set_trace()

        #get consensus interaction
        consensusData = self.calcEnergySequence(self.consensus,"consensus", mutPair)
        consensusInteraction = self.getInteraction(consensusData[4],consensusData[5],consensusData[3])
        print(consensusData, consensusInteraction)

        ##mutationpair, seq index, delta delta, d12, d1, d2

        index = 0
        for seqArr in valueStack:
            seq = self.seqStack[seqArr[1]]
            actualDeltam1m2 = seqArr[3]
            Dm1 = seqArr[4]
            Dm2 = seqArr[5]
            DDe = seqArr[2]
            hd = self.getHDNum(seq)
            seqMutList =  len(hd)

            comparisonData = self.getComparisonDelta(Dm1,Dm2,actualDeltam1m2)

            distFromCompValue = comparisonData[0]
            isRescue = comparisonData[1]

            ##check interaction
            ##seqData = [mutPair,index,interaction,distFromMaxDelta1,actualDeltam1m2]
            seqData = [Dm1,Dm2,actualDeltam1m2,DDe,hd,seqMutList,distFromCompValue]


            ##new comparison func
            index+=1
            if isRescue:
                rescue.append(seqData)
            elif distFromCompValue >= 0:
                compensate.append(seqData)
            else:
                antag.append(seqData)
        
        return [rescue,compensate,antag]


    ##gets comparison data. if rescue, gets max distance from top mutation. Else, gets distance from bottom mutation.
    def getComparisonDelta(self,Dm1,Dm2,Dm1m2):
        compVal1 = max(Dm1, Dm2)
        compVal2 = min(Dm1, Dm2)
        isRescue = False
        if Dm1m2 > compVal1:
            delta = Dm1m2 - compVal1
            isRescue = True
        else:
            delta = Dm1m2 - compVal2
        return [delta,isRescue]
    
    def importValueStacks(self,filename):
        with open(filename) as oFile:
            returnArr = json.load(oFile)

        self.importedStacks.append(returnArr)        


    #returns interaction
    def getInteraction(self,Dm1,Dm2, actual):
        projectedDelta = Dm1 + Dm2
        if actual < projectedDelta:
            return self.INTERACTION3
        elif actual > max(Dm1,Dm2):
            return self.INTERACTION1
        else:
            return self.INTERACTION2
   
        


    ##compares value stack interactions with consensus. sorts into one big list by interaction, delta difference.
    ##combine all relevant data into one entry. 
    ## label,index,interaction,consensusInteraction,deltaDifference
    ##returns 2 arrays:[same as consensus, diff to consensus]
    def compareToConsensus(self,loadValueStacks):
        ##gets values: [[mutation,computed stack values],[mutation,computed stack values],...]
        sequenceCalculations = []
        for valueStack in self.importedStacks:
            sequenceCalculations += valueStack
        
            ##sequenceCalulations: [[mutation,computed stack values],[mutation,computed stack values]]
        for mutationData in sequenceCalculations:
            mutName = self.arrToString(mutationData[0][0]) + "-" + self.arrToString(mutationData[0][1])
            print(mutName)
            seqData = self.getInteractionsAndDifference(mutationData)
            rescue = seqData[0]
            compensate = seqData[1]
            antag = seqData[2]
            

            #sorts
            rescue = sorted(rescue,key=lambda x: x[6],reverse=True)
            compensate = sorted(compensate,key=lambda x: x[6],reverse=True)
            antag = sorted(antag,key=lambda x: x[6],reverse=True)

            print("computed for mutation ",mutName)
            print("sorted by distance from Dm1m2 to highest-fitness single mutation")
            ##             seqData = [Dm1,Dm2,actualDeltam1m2,DDe,hd,seqMutList,distFromCompValue]


            # [mutationArray,name,deltaDeltaE,m1m2Delta,m1Delta,m2Delta]
            consensusRow = self.calcEnergySequence(self.consensus,"consensus", mutationData[0])

            lastRow = [consensusRow[4],consensusRow[5],consensusRow[3],consensusRow[2],"None",0,0]

            rescue.insert(0,lastRow)
            compensate.insert(0,lastRow)
            antag.insert(0,lastRow)

            #writes files
            with open("out/"+mutName+".rescue.json", "w") as outfile1:
                json.dump(rescue, outfile1)
            with open("out/"+mutName+".compensate.json", "w") as outfile2:
                json.dump(compensate, outfile2)
            with open("out/"+mutName+".antag.json", "w") as outfile3:
                json.dump(antag, outfile3)


            self.writeToCSV(rescue,mutName+"rescue")
            self.writeToCSV(compensate,mutName+"compensate")
            self.writeToCSV(antag,mutName+"antag")

    def writeToCSV(self,data,name):
        fields = ["Dm1","Dm2","Dm1m2","DDe","list of mutations","HD","distance from max/min D single mutation"]
        with open("out/"+name+".csv", "w") as outfile1:                
            write = csv.writer(outfile1)            
            write.writerow(fields)
            write.writerows(data)

    def createValueStacks(self):
        with open("out/newValueStack.json","w") as outfile:
            json.dump(self.computeValueStacks(),outfile)
        
        
        
        

            


    
    


        










        

        