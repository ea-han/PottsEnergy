import json
import os


##this class imports sequences+mutation pair data from CompileSequences and calculates their energies.
class CreateValueStacks:

    ## [mutation pair, list of sequences]
    importedStacks = []
    jFileArray = []
    reduxName = ''

    def __init__(self,jFileArray,reduxName):
        self.importedStacks = []
        self.jFileArray = jFileArray
        self.reduxName = reduxName

    #loads all stacks fro src/compiledSeqs
    def loadStacksFromJson(self):
        for filename in os.scandir("src/compiledSeqs"):
            if filename.is_file():

                print("found file: ", filename.path)
                with open(filename.path) as iFile:
                     self.importedStacks.append(json.load(iFile))   



    def arrToString(self,mutation):
        altMutations = ','.join(mutation[2])
        newString = ''
        newString = newString + str(mutation[1])
        newString = newString + str(mutation[0]+1)
        newString = newString + altMutations
        return newString

    def mutPairToString(self,mutationPair):
        string1 = self.arrToString(mutationPair[0])
        string2 = self.arrToString(mutationPair[1])
        return string1 + "-" + string2


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
    
    ##calculates for sequence by reversing the mutation etc.
    ##NEW
    def calcEnergySequenceReversed(self, sequence,name,mutationArray): 

        defaultEnergy = self.getEnergy(self.mutateSequence2Reverse(sequence,mutationArray[0],mutationArray[1])) ##works
        m2Energy = self.getEnergy(self.mutateSequence1Reverse(sequence, mutationArray[0]))
        m1Energy = self.getEnergy(self.mutateSequence1Reverse(sequence, mutationArray[1]))
        m1m2Energy = self.getEnergy(sequence) ##works

        m1Delta = defaultEnergy - m1Energy 
        m2Delta = defaultEnergy - m2Energy
        m1m2Delta = defaultEnergy - m1m2Energy ##works
        deltaDeltaE = m1m2Delta - (m1Delta + m2Delta)

        m1Delta += deltaDeltaE
        m2Delta += deltaDeltaE

        ##name is index
        return [mutationArray,name,deltaDeltaE,m1m2Delta,m1Delta,m2Delta]
    

        ##undoes the mutation
    def mutateSequence1Reverse(self, seq,mutation):
        sequence = list(seq)
        sequence[mutation[0]] = mutation[1]
        return sequence

    ##undoes both mutations
    def mutateSequence2Reverse(self, seq,mutation1,mutation2):
        sequence = list(seq)
        sequence[mutation1[0]] = mutation1[1]
        sequence[mutation2[0]] = mutation2[1]
        return sequence

    # This method calculates value stacks for all mutations and outputs them to json.
    # imports categorized sequences per mutation and calculates energies.
    # NEW
    def createValueStacks(self):
        self.loadStacksFromJson()

        outputValueStack = []

        for data in self.importedStacks:
            mutPair = data[0]
            mutPairReduced = self.reduceMutations(mutPair)

            reducedStack = data[1]
            energiesStack = []


            for sequence in reducedStack:
                energiesStack.append(self.calcEnergySequenceReversed(sequence,self.mutPairToString(mutPair),mutPairReduced))

            print("appended 1 value stack")
            stack1 = [mutPair, energiesStack]
            self.writeToJson(stack1,"src/stackjson/"+self.mutPairToString(mutPair)+".stack")

            outputValueStack.append(stack1)

        self.writeToJson(outputValueStack,"src/allValueStacks")
                

    def writeToJson(self,data,name):
        with open(name+".json", "w") as outfile:
                json.dump(data, outfile)

    def reduceMutations(self,mutationPair): 
        ctr = 0
        mutation1 = mutationPair[0]
        mutation2 = mutationPair[1]
        with open(self.reduxName) as kFile:
            for line in kFile:
                if ctr == mutation1[0]-1:
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
                    
                    
                    
                elif ctr == mutation2[0]-1:
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
