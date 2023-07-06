import numpy as np
import io as io
import sys
import EnergyAverage

def create(index1, letter1, letterarr1, index2, letter2, letterarr2):
    index1 -= 1
    index2 -= 1

    ##start
    jFile = open("in.fullseq")
    totalSeq = []
    lines = jFile.readlines()
    for line in lines:
        strHolder = line.strip()
        appendList = list(strHolder)
        totalSeq.append(appendList)

    cFile = open("in.consensus.reduce4.seq", 'r')
    consensus = cFile.readline().strip()
    
    mutation1 = [index1,letter1,letterarr1]
    mutation2 = [index2,letter2,letterarr2]

    jFileArray = np.load("J.npy")
    reduxName = "in.reduce4.redux"

    energy1 = EnergyAverage.EnergyAverage(totalSeq, consensus, mutation1, mutation2, jFileArray, reduxName)
    energy1.calcAverage()

def main1():
    inputArray =[[140,'G',['S'],148,'Q',['H']],
                 [6,'D',['G','S','L','T'],7,'K',['E']],
                 [143,'Y',['C','W','T'],230,'S',['E','R']],
                 [157,'E',['F','G','K'],160,'K',['A','G','M','Q']],
                 [119,'S',['V','E','G'],122,'T',['C','I']],
                 [138,'E',['K','T'],148,'Q',['N']],
                 [219,'K',['N','H'],222,'N',['C','K']],
                 [110,'V',['A','E','I'],142,'P',['A']],
                 [32,'V',['C','I'],222,'N',['A','F','H']],
                 [227,'Y',['H','F','A'],253,'D',['A','N','F','Y']],
                 [227,'Y',['H','F','A'],253,'D',['C','H']],
                 [70,'G',['A','S'],119,'S',['A','D','H']],
                 [157,'E',['F','K','G'],160,'K',['I','V','P']],
                 [60,'I',['C','E','G'],182,'I',['C','E','G']]
                 ]

    for input in inputArray:
        create(input[0],input[1],input[2],input[3],input[4],input[5])

main1()  
