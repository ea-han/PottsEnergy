import Energy
import io

def main():
    seqFile = input("input sequence filename: \n")
    jFile = input("input npy filename\n")
    eOb = Energy.Energy(seqFile, jFile)
    eOb.simplify()
    eOb.getEnergy()
    print("energy = ", eOb.energy)

main()
