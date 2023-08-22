import numpy as np
import io as io
import sys
import EnergyAverage2
import pandas as pd
import pdb

sheetArray = []
plotArrayMut = []
plotArrayC = []
plotArray1 = []
plotArray2 = []
plotArray3 = []


##seqStack, consensus, mutationList, jFileArray, reduxName
def create(seqStackfile, mutationList):
    ##load files
    jFile = open(seqStackfile)
    seqStack = []
    lines = jFile.readlines()
    for line in lines:
        strHolder = line.strip()
        appendList = list(strHolder)
        seqStack.append(appendList)

    cFile = open("src/IN_CONSENSUS", 'r')
    consensus = cFile.readline().strip()

    jFileArray = np.load("src/J.npy")
    reduxName = "src/in.reduce4.redux"

    ##subtract one from every index in mutationlist
    for mutationPair in mutationList:
        mutationPair[0][0] = mutationPair[0][0] - 1
        mutationPair[1][0] = mutationPair[1][0] - 1

    ##init energy2
    energy = EnergyAverage2.EnergyAverage2(seqStack, consensus, mutationList, jFileArray, reduxName)
    return energy

##Runsc
def loadFromFile(fileName):
    ##INPUT FORMAT: src,idx,dest1,dest2,...,destn
    ##              src,idx,dest1,dest2,...,destn
    ##              \n
    inputArr = []
    insertArr1 = []
    insertArr2 = []


    file = open(fileName)
    lines = file.readlines()
    line1 = None
    line2 = None
    for line in lines:
        if line == "\n":
            ##add to the array
            line1 = line1.strip()
            line1 = line1.split(',')
            insertArr1.append(line1.pop(1))
            insertArr1.append(line1.pop(0))
            insertArr1.append(line1)

            line2 = line2.strip()
            line2 = line2.split(',')
            insertArr2.append(line2.pop(1))
            insertArr2.append(line2.pop(0))
            insertArr2.append(line2)

            insertArr1[0] = int(insertArr1[0])
            insertArr2[0] = int(insertArr2[0])

            inputArr.append([insertArr1,insertArr2])
            print("appended ",insertArr1,insertArr2)
            insertArr1 = []
            insertArr2 = []

            line1 = None
            line2 = None 

        elif line1 == None:
            line1 = line

        elif line2 == None:
            line2 = line
        
    return inputArr



def main():
    inputArray = loadFromFile("src/testinput")
    object1 = create("src/in.fullseq",inputArray)
    ##object1.preReduceToJson()
    object1.loadJson()
    object1.importValueStacks("src/valueStack.json")
    object1.reduceAllMutations()
    object1.compareToConsensus(True)

    

main()




