import unittest
import EnergyAverage
import numpy as np

class TestEnergy(unittest.TestCase):
    ##reducednl43
    nl43Reduced = "ABAAABCBACABAACBDBABBDBDBBBAAACAABAAABDAABAABAAABCBAAAAACACBAABAADACCBADAACACABBABABAAABABBAABDAABCAAAABABAAAABDDBCAAADABCCCDDAAAAAAAABDBDACBADAAABDABBAAAACBBADAAACCAAABAAABAAABAABABAABABABAAAACCAAAABBAABBBBCACDCAACDBDCCCAAAAABADDBCAAAAAAABAABAAAAAAABBDBDBACCBBBA"
    inputStringNL43 = list("FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTVHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVIESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRDPVWKGPAKLLWKGEGAVVIQDNSDIKVVPRR")
    nl43Reduced = list(nl43Reduced)

    def setUp(self):
        input =    [143,'Y',['A'],148,'Q',['Q']]
        return self.create(input[0],input[1],input[2],input[3],input[4],input[5],"in.nl43.seq")


    def setUpAlternate(self):
        input =    [140,'G',['G'],143,'Y',['A']]
        return self.create(input[0],input[1],input[2],input[3],input[4],input[5],"in.nl43.seq")

    def create(self,index1, letter1, letterarr1, index2, letter2, letterarr2, seqfile):
        index1 -= 1
        index2 -= 1
        ##start
        jFile = open(seqfile)
        totalSeq = []
        lines = jFile.readlines()
        for line in lines:
            strHolder = line.strip()
            appendList = list(strHolder)
            totalSeq.append(appendList)

        cFile = open("IN_CONSENSUS", 'r')
        consensus = cFile.readline().strip()
        
        mutation1 = [index1,letter1,letterarr1]
        mutation2 = [index2,letter2,letterarr2]

        jFileArray = np.load("J.npy")
        reduxName = "in.reduce4.redux"

        energy1 = EnergyAverage.EnergyAverage(totalSeq, consensus, mutation1, mutation2, jFileArray, reduxName)
        return energy1


    def testReduce(self):##passed
        newOb = self.setUp()
        self.assertEqual(newOb.reduceSequence(self.inputStringNL43), self.nl43Reduced)
    
    # Tests NL4-3 Q148H single mutation with calconetooutput
    def testCalcOneToOutputReference(self):
        newOb = self.setUp()
        nl43Arr = newOb.calcOneToOutput(self.inputStringNL43,"nl43")
        dde = nl43Arr[2]
        dm1m2 = nl43Arr[3]
        dm1 = nl43Arr[4]
        dm2 = nl43Arr[5]
        print(nl43Arr)
        self.assertEquals(dm1m2,-9.916907684877514839e+00)
        self.assertEquals(dde,dm1m2-(dm1+dm2))

        newOb = self.setUpAlternate()
        nl43Arr = newOb.calcOneToOutput(self.inputStringNL43,"nl43")
        dde = nl43Arr[2]
        dm1m2 = nl43Arr[3]
        dm1 = nl43Arr[4]
        dm2 = nl43Arr[5]
        print(nl43Arr)
        self.assertEquals(dm1m2,-6.705908505246043205e+00)
        self.assertEquals(dde,dm1m2-(dm1+dm2))

        ##"dm1m2" is actually Delta M1 because Dm2 is 0 
        

