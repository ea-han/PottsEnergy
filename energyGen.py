import Energy
import io
def gen():
    eOb = Energy.Energy("IN.seq", "J.npy")
    eOb.simplify()
    eOb.genNpy(4,"4.npy")

gen()