import Energy

def average():

    energyOb = Energy.Energy("IN_CONSENSUS","J.npy","input.txt","in.reduce4.redux","outputtest_fullseq.txt")
    ##nFile = open("in.fullseq", 'r')
    nFile = open("in.fullseq", 'r')

    inputList = []
    lines = nFile.readlines()
    ctr = 0

    ##reads line by line and adds sequence array to inputlist
    for line in lines:
        strHolder = line.strip()
        appendList = list(strHolder)
        inputList.append([])
        inputList[ctr] = appendList
        ctr += 1
    energyOb.mutateAverageOf3(inputList)



def consensus():
    energyOb = Energy.Energy("IN_CONSENSUS","J.npy","input.txt","in.reduce4.redux","outputtest_consensus.txt")
    energyOb.mutate()


def main():
    done = False
    OUTPUT_NAME = "0output.txt"
    seqFile = input("input sequence filename: \n")
    outputName = OUTPUT_NAME
    ctr = 1
    while not done:
        jFile = input("input npy filename\n")
        mutList = input("input mutation list:\n")
        reduxName = input("input .redux filename:")
        energyOb = Energy.Energy(seqFile,jFile,mutList,reduxName,outputName)
        energyOb.mutate()
        outputName = outputName[1:]
        outputName = str(ctr) + outputName
        seqFile = input("Input new sequence filename or q to quit:\n") 
        if seqFile == "q":
            done = True
        ctr+=1


##Input: wildtype sequence, list of input sequences
##get the positions of the mutation, find the mutations in a set of sequences, calculate delta E for m1m2, m1, m2. subtract m1+m2 from m1m2, get the average for the list of sequencs

##First, get inputs and calculate delta-delta e

##main()
average()
consensus()