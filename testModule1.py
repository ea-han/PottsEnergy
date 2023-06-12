import unittest
import Energy

class TestEnergy(unittest.TestCase):
    energyOb = None
    unreducedString = "ABAAABCBACABAACBDBABBDBDBBBAAACAABAAABDAABAABAAABCBAAAAACACBAABAADACCBAAAACACABBABABAAABABBAABDAABCAAAABABAAAABDCBCAAADABCCCDDAAAAAAAABDBDACBADAAABDABAAAAACBBADAAACCAAABAAABAAABAABABAABABABAAAACCAAAABBAABBBBCACDCAACDBDCCCAAAAABADDBCACAAAAABAABAAAAAAABBDBDBACCBBBA"

    wildType = list(unreducedString)
    def setUp(self):
        self.energyOb = Energy.Energy("unreduced.seq","J.npy","test.txt","in.reduce4.redux","testOutput.txt")

    def testReduce(self):
        self.assertEqual(self.energyOb.reduce(self.energyOb.nVals,"in.reduce4.redux"), self.wildType)