import unittest
import Energy

class TestEnergy(unittest.TestCase):
    energyOb = None
    unreducedString = "ABAAABCBACABAACBDBABBDBDBBBAAACAABAAABDAABAABAAABCBAAAAACACBAABAADACCBADAACACABBABABAAABABBAABDAABCAAAABABAAAABDDBCAAADABCCCDDAAAAAAAABDBDACBADAAABDABBAAAACBBADAAACCAAABAAABAAABAABABAABABABAAAACCAAAABBAABBBBCACDCAACDBDCCCAAAAABADDBCAAAAAAABAABAAAAAAABBDBDBACCBBBA"

    wildType = list(unreducedString)
    def setUp(self):
        self.energyOb = Energy.Energy("in.nl43.seq","J.npy","test.txt","in.nl43.reduce4.redux","testOutput.txt")

    def testReduce(self):
        self.assertEqual(self.energyOb.reduce(self.energyOb.nVals,"in.reduce4.redux"), self.wildType)