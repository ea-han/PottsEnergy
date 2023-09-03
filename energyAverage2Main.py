import numpy as np
import io as io
import sys
import pandas as pd
import json
import ClassifySequences
import CreateValueStacks
import CompileSequences

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
    energy = ClassifySequences.ClassifySequences(seqStack, consensus, mutationList, jFileArray, reduxName)
    return energy

def makeDataStacks(seqStackFile,mutationList):
    jFile = open(seqStackFile)
    seqStack = []
    lines = jFile.readlines()
    for line in lines:
        strHolder = line.strip()
        appendList = list(strHolder)
        seqStack.append(appendList)

    ds = CompileSequences.CompileSequences("src/in.reduce4.redux",seqStack,mutationList)
    ds.computeAllStackPairs()


def makeSeqStacks():
    jArray = np.load("src/J.npy")
    ob = CreateValueStacks.CreateValueStacks(jArray,"src/in.reduce4.redux")
    ob.createValueStacks()


##Runsc
def loadFromFile(fileName):
    ##INPUT FORMAT: src,idx,dest1,dest2,...,destn
    ##              src,idx,dest1,dest2,...,destn
    ##              \n
    inputArr = []
    insertArr1 = []
    insertArr2 = []


    file = open(fileName,"r")
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

def foo():
    inputArray = loadFromFile("src/testinput")
    object1 = create("src/in.fullseq",inputArray)
    print(''.join(object1.consensus))

def main():
    inputArray = loadFromFile("src/testinput")
    ##object1 = makeDataStacks("src/in.fullseq",inputArray)
    ##makeSeqStacks()
    ob2 = create("src/in.fullseq",inputArray)
    ob2.compareToConsensus()
    
##foo()
main()





