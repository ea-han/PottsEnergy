import json
## This class searches through seqStack for sequences that contain the 
## desired mutation pairs and outputs categorized,reduced sequences.
class CompileSequences:
    reduxName = ""
    seqStack = []
    mutList = []

        ## [mutation pair, list of sequences]
    def __init__(self,reduxName,seqStack,mutList):
        self.reduxName = reduxName
        self.seqStack = seqStack
        self.mutList = mutList



    def arrToString(self,mutation):
        altMutations = ','.join(mutation[2])
        newString = ''
        newString = newString + str(mutation[1])
        newString = newString + str(mutation[0])
        newString = newString + altMutations
        return newString

    def mutPairToString(self,mutationPair):
        string1 = self.arrToString(mutationPair[0])
        string2 = self.arrToString(mutationPair[1])
        return string1 + "-" + string2


    def writeToJson(self,data,name):
        with open(name+".json", "w") as outfile:
                json.dump(data, outfile)

    ##searches through seqStack for single mutation pair
    def searchForMutPair(self,mutPair):
        print(mutPair)
        mutation1 = mutPair[0]
        mutation2 = mutPair[1]
        index1 = mutation1[0]-1
        index2 = mutation2[0]-1

        

        foundSeqs = []
        ctr = 0

        for seq in self.seqStack:
            if (seq[index1] in mutation1[2]) and (seq[index2] in mutation2[2]):
                foundSeqs.append(seq)
            ctr+=1

        if (len(foundSeqs) == 0):
            print("no sequences found for ",mutPair)
            return -1      
        
        return [mutPair,foundSeqs]

    ##computes,reduces and exports a stackpair for a single mutationPair
    def getStackPair(self,mutPair):
        stackPair = self.searchForMutPair(mutPair)

        if stackPair == -1:
            return -1
        
        newStack = []
        stack = stackPair[1]
        for seq in stack:
            seq = self.reduceSequence(seq)
            newStack.append(seq)

        newStackPair = [mutPair,newStack]

        self.writeToJson(newStackPair,"src/compiledSeqs/"+self.mutPairToString(mutPair))
        return newStackPair

    ##this method computes all stackpairs.
    def computeAllStackPairs(self):
        fullStack = []
        for mutPair in self.mutList:
            stackPair = self.getStackPair(mutPair)
            if (stackPair != -1):
                fullStack.append(stackPair)

        print("compiled all stacks")
        self.writeToJson(fullStack,'src/allCompiledSeqs')




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
