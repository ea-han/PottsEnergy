import re
def strip(inputString):
    inputString = inputString.replace('\n','')
    inputString = inputString.replace(' ','')
    inputString = inputString.upper()
    nums = r'[0-9]'
    inputString = re.sub(nums,'',inputString)
    inputString = inputString[0:263]
    return inputString

def main():
    file = open("stringInput.txt")
    output = open("outputSeqs.txt",'w')
    lines = file.readlines()
    inputString = ""
    for line in lines:
        if line != "\n":
            inputString += line
        else:
            addString = strip(inputString)
            output.write(addString)
            output.write('\n')
            inputString = ""




main()