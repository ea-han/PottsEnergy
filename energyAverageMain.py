import numpy as np
import io as io
import sys
import EnergyAverage
import matplotlib as mpl
import pandas as pd
mpl.use('TkAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.random import rand
from pylab import figure
import plotly.express as px
import plotly.graph_objects as go

sheetArray = []
plotArrayMut = []
plotArrayC = []
plotArray1 = []
plotArray2 = []
plotArray3 = []

mutArr1 = [
"G140S-Q148H",
"Y143C/W/T-S230E/R",
"G140S-Q148R",
"E138K/T-Q148N/K",
"V110A/E/I-P142A",
"G70A/S-S119A/D/H",
"T97A-Y143R"
]
mutArr2 = [
"G140S-Q148H",
"D6G/S/L/T-K7E",
"Y143C/W/T-S230E/R",
"E157F/G/K-K160A/G/M/Q",
"S119V/E/G-T122C/I",
"G140S-Q148R",
"E138K/T-Q148N/K",
"K219N/H-N222C/K"
]

def create(index1, letter1, letterarr1, index2, letter2, letterarr2, seqfile):
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

    cFile = open("src\IN_CONSENSUS", 'r')
    consensus = cFile.readline().strip()
    
    mutation1 = [index1,letter1,letterarr1]
    mutation2 = [index2,letter2,letterarr2]

    jFileArray = np.load("src/J.npy")
    reduxName = "src\in.reduce4.redux"

    energy1 = EnergyAverage.EnergyAverage(totalSeq, consensus, mutation1, mutation2, jFileArray, reduxName)
    return energy1

def loadFromFile(fileName):
    ##INPUT FORMAT: src,idx,dest1,dest2,...,destn
    ##              src,idx,dest1,dest2,...,destn
    ##              \n
    inputArr = []
    insertArr = []

    file = open(fileName)
    lines = file.readlines()
    line1 = None
    line2 = None
    for line in lines:
        if line == "\n":
            ##add to the array
            line1 = line1.strip()
            line1 = line1.split(',')
            insertArr.append(line1.pop(1))
            insertArr.append(line1.pop(0))
            insertArr.append(line1)

            line2 = line2.strip()
            line2 = line2.split(',')
            insertArr.append(line2.pop(1))
            insertArr.append(line2.pop(0))
            insertArr.append(line2)

            insertArr[0] = int(insertArr[0])
            insertArr[3] = int(insertArr[3])

            inputArr.append(insertArr)
            print(insertArr)
            insertArr = []
            line1 = None
            line2 = None 

        elif line1 == None:
            line1 = line

        elif line2 == None:
            line2 = line
        
    print(inputArr)
    return inputArr


def buildPlotly2d():
    thdArray = []
    thdArray.append(plotArrayC)
    thdArray.append(plotArray1)
    thdArray.append(plotArray2)    
    thdArray.append(plotArray3)

    nameArr = ["Consensus","HXB2","LAI-IIIB","NL4-3"]
    i = 0

    for plotArr in thdArray:
        inputArr = []
        for elem in plotArr:
            inputArr.append([elem[0],"DeltaM1M2",elem[3],elem[2]])
            inputArr.append([elem[0],"DeltaM1",elem[4],elem[2]])
            inputArr.append([elem[0],"DeltaM2",elem[5],elem[2]])
            inputArr.append([elem[0],"ProjectedM1M2",(elem[4]+elem[5]),elem[2]])
        plotly2d(inputArr,nameArr[i])
        print("and ",plotArr[0][1])
        i+=1

def mainAll():
    inputStringHXB2 = list("FLDGIDKAQDEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTSATVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRR")
    inputStringLAIIIIB = list("FLDGIDKAQDEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTSATVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRNPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRR")
    inputStringNL43 = list("FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTVHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVIESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRDPVWKGPAKLLWKGEGAVVIQDNSDIKVVPRR")
    

    inputArray = loadFromFile("src/testinput")





    for input in inputArray:
            energy1 = create(input[0],input[1],input[2],input[3],input[4],input[5],"HXB2")
            energy2 = create(input[0],input[1],input[2],input[3],input[4],input[5],"LAI_IIIB")
            energy3 = create(input[0],input[1],input[2],input[3],input[4],input[5],"NL4-3")
            HXB2Arr = energy1.calcOneToOutput(inputStringHXB2,"HXB2")
            LAIIIBArr = energy2.calcOneToOutput(inputStringLAIIIIB,"LAI-IIIB")
            NL43Arr = energy3.calcOneToOutput(inputStringNL43,"NL4-3")
            ConArr = energy1.calcConsensus()
            ConArr.insert(1,"Consensus")

            plotArrayMut.append(ConArr[0])
            plotArrayC.append(ConArr)
            plotArray1.append(HXB2Arr)
            plotArray2.append(LAIIIBArr)
            plotArray3.append(NL43Arr)

            sheetArray.append(ConArr)
            sheetArray.append(HXB2Arr)
            sheetArray.append(LAIIIBArr)
            sheetArray.append(NL43Arr)



def plotly3d():
    df2 = pd.DataFrame(sheetArray,index=None,columns=['Mutation Pair','Strain','DeltaDelta','DeltaM1M2','DeltaM1','DeltaM2'])
    df2.drop(0)
    df2.drop(1)
        
    fig = px.scatter_3d(df2, x="DeltaM1", y="DeltaM2",z="DeltaM1M2", color="Strain",
                 labels={
                     "DeltaM1": "Delta M1",
                     "DeltaM2": "Delta M2",
                     "DeltaM1M2": "Delta M1-M2",
                     "Strain": "Strain"
                 },
                title=":)")

    fig = figure()
    ax = fig.add_subplot(projection='3d')

    for i in range(len(plotArrayC)): #plot each point + it's index as text above
        ax.scatter(plotArrayC[i][4],plotArrayC[i][5],plotArrayC[i][3],label="Consensus",color='b') 
        if labels:
            ax.text(plotArrayC[i][4],plotArrayC[i][5],plotArrayC[i][3],  '%s' % (plotArrayC[i][0] + ' ' + plotArrayC[i][1]), size=5, zorder=1,  
        color='k') 

    for i in range(len(plotArray1)): #plot each point + it's index as text above
        ax.scatter(plotArray1[i][4],plotArray1[i][5],plotArray1[i][3],color='r',label="HXB2",marker='o') 
        if labels:
            ax.text(plotArray1[i][4],plotArray1[i][5],plotArray1[i][3],  '%s' % (plotArray1[i][0] + ' ' + plotArray1[i][1]), size=5, zorder=1,  
        color='k') 

    for i in range(len(plotArray2)): #plot each point + it's index as text above
        ax.scatter(plotArray2[i][4],plotArray2[i][5],plotArray2[i][3],color='g',label="LAI III-B",marker='^') 
        if labels:
            ax.text(plotArray2[i][4],plotArray2[i][5],plotArray2[i][3],  '%s' % (plotArray2[i][0] + ' ' + plotArray2[i][1]), size=5, zorder=1,  
        color='k') 

    for i in range(len(plotArray3)): #plot each point + it's index as text above
        ax.scatter(plotArray3[i][4],plotArray3[i][5],plotArray3[i][3],color='y',label='NL4-3',marker='*') 
        if labels:
            ax.text(plotArray3[i][4],plotArray3[i][5],plotArray3[i][3],  '%s' % (plotArray3[i][0] + ' ' + plotArray3[i][1]), size=5, zorder=1,  
        color='k') 

    ax.set_xlabel('Delta M1')
    ax.set_ylabel('Delta M2')
    ax.set_zlabel('Delta M1-M2')

    ax.set_xlim([-15,5])    
    ax.set_ylim([-15,5])    
    ax.set_zlim([-15,5])    

    

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())



    plt.show()

def plotly2d(inputArr,name):
    df2 = pd.DataFrame(inputArr,index=None,columns=["MutationPair","DeltaType","Energy","DDE"])
    df2['symbol'] = pd.Series(["1" for x in range(len(df2.index))])
    print(df2)


    x_end = []
    y_end = []
    x_start = []
    y_start = []    

    for idx in df2.index:
        if df2["DeltaType"][idx] == "ProjectedM1M2":
            x_start.append(df2["MutationPair"][idx])
            y_start.append(df2["Energy"][idx])
            df2.iat[idx,4] = "2"

        if df2["DeltaType"][idx] == "DeltaM1M2":
            x_end.append(df2["MutationPair"][idx])
            y_end.append(df2["Energy"][idx])
            df2.iat[idx,4] = "3"


    df2 = df2.sort_values(['DDE'])
    

    fig = px.scatter(df2, symbol="symbol",y="Energy", x="MutationPair", color="DeltaType", 
                 labels={
                     "Energy": "Delta E",
                     "MutationPair": "Mutation Pair",
                     "DeltaType": "Mutations"
                 },
                title=name)
    fig.update_layout(yaxis_range=[-25,5])
    fig.update_traces(marker=dict(size=12))

    print(df2)

    list_of_all_arrows = []
    for x0,y0,x1,y1 in zip(x_end, y_end, x_start, y_start):
        arrow = go.layout.Annotation(dict(
                        x=x0,
                        y=y0,
                        xref="x", yref="y",
                        text="",
                        showarrow=True,
                        axref="x", ayref='y',
                        ax=x1,
                        ay=y1,
                        arrowhead=3,
                        arrowwidth=1.5,
                        arrowcolor='purple',)
                    )
        list_of_all_arrows.append(arrow)

    fig.update_layout(annotations=list_of_all_arrows)
    print(df2)
    fig.write_html(name + ".html")
    fig.show()
    print("sanity check: name = ",name)


def plot2d(mutArr):
    import matplotlib.pyplot as plt
    import numpy as np; np.random.seed(0)
    import matplotlib.transforms as transforms
    
    SMALL_SIZE = 7
    MEDIUM_SIZE = 12
    BIGGER_SIZE = 14

    plt.yticks(fontsize=7)
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=6)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    plotArrayCCopy = []
    plotArray1Copy = []
    plotArray2Copy = []
    plotArray3Copy = []

    for elem in plotArrayC:
        print(elem[0])
        print(mutArr)
        if elem[0] in mutArr:
            plotArrayCCopy.append(elem)
            
    for elem in plotArray1:
        if elem[0] in mutArr:
            plotArray1Copy.append(elem)
            
    for elem in plotArray2:
        if elem[0] in mutArr:
            plotArray2Copy.append(elem)
            
    for elem in plotArray3:
        if elem[0] in mutArr:
            plotArray3Copy.append(elem)


    values = np.random.rand(300, 3)

    plt.figure()

    offset = lambda p: transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
    trans = plt.gca().transData

    for i in range (0,len(plotArrayCCopy)):
        sc1 = plt.scatter(plotArrayCCopy[i][0], plotArrayCCopy[i][4], c = 'red', s = 25, transform=trans+offset(-10),label="Consensus M1",marker="x")##Consensus
        plt.scatter(plotArray1Copy[i][0], plotArray1Copy[i][4], c = 'blue', s = 25,marker="x",label="HXB2 M1")
        plt.scatter(plotArray2Copy[i][0], plotArray2Copy[i][4], c = 'green', s = 25, transform=trans+offset(10),marker="x",label="LAI-IIIB M1")
        plt.scatter(plotArray3Copy[i][0], plotArray3Copy[i][4], c = 'orange', s = 25, transform=trans+offset(20),marker="x",label="NL4-3 M1")
        
        plt.scatter(plotArrayCCopy[i][0], plotArrayCCopy[i][5], c = 'red', s = 25, transform=trans+offset(-10),marker="*",label="Consensus M2")##Consensus
        plt.scatter(plotArray1Copy[i][0], plotArray1Copy[i][5], c = 'blue', s = 25,label="HXB2 M2",marker="*")
        plt.scatter(plotArray2Copy[i][0], plotArray2Copy[i][5], c = 'green', s = 25, transform=trans+offset(10),marker="*",label="LAI-IIIB M2")
        plt.scatter(plotArray3Copy[i][0], plotArray3Copy[i][5], c = 'orange', s = 25, transform=trans+offset(20),marker="*",label="NL4-3 M2")

        plt.scatter(plotArrayCCopy[i][0], plotArrayCCopy[i][3], c = 'red', s = 25, transform=trans+offset(-10),marker="^",label="Consensus M1-M2")##Consensus
        plt.scatter(plotArray1Copy[i][0], plotArray1Copy[i][3], c = 'blue', s = 25,marker="^",label="HXB2 M1-M2")
        plt.scatter(plotArray2Copy[i][0], plotArray2Copy[i][3], c = 'green', s = 25, transform=trans+offset(10),marker="^",label="LAI-IIIB M1-M2")
        plt.scatter(plotArray3Copy[i][0], plotArray3Copy[i][3], c = 'orange', s = 25, transform=trans+offset(20),label="NL4-3 M1-M2",marker="^")

        plt.scatter(plotArrayCCopy[i][0], plotArrayCCopy[i][4] + plotArrayCCopy[i][5], c = 'red', s = 25, transform=trans+offset(-10),marker="o",label="Projected Consensus M1-M2")##Consensus
        plt.scatter(plotArray1Copy[i][0], plotArray1Copy[i][4] + plotArray1Copy[i][5], c = 'blue', s = 25,marker="o",label="Projected HXB2 M1-M2")
        plt.scatter(plotArray2Copy[i][0], plotArray2Copy[i][4] + plotArray2Copy[i][5], c = 'green', s = 25, transform=trans+offset(10),label="Projected LAI-IIIB M1-M2",marker="o")
        plt.scatter(plotArray3Copy[i][0], plotArray3Copy[i][4] + plotArray3Copy[i][5], c = 'orange', s = 25, transform=trans+offset(20),marker="o",label="Projected NL4-3 M1-M2")

        plt.text(i-.1,2,"- Consensus",rotation = 90,
         rotation_mode = 'anchor',
         transform_rotates_text = True)
        plt.text(i,2,"- HXB2",rotation = 90,
         rotation_mode = 'anchor',
         transform_rotates_text = True)
        plt.text(i+.1,2,"- LAI-IIIB",rotation = 90,
         rotation_mode = 'anchor',
         transform_rotates_text = True)
        plt.text(i+.2,2,"- NL4-3",rotation = 90,
         rotation_mode = 'anchor',
         transform_rotates_text = True)

        
        
    plt.xlabel("Mutations")
    plt.ylabel("Delta E")
    plt.ylim([-20,8])


    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))


    plt.xticks(fontsize=5,rotation=60)
    plt.tick_params('both',length=5,width=2,labelsize=5)
    plt.legend(by_label.values(), by_label.keys(),bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)

    plt.show()


mainAll()
df = pd.DataFrame(sheetArray,index=None,columns=['Mutation Pair','','DeltaDelta','DeltaM1M2','DeltaM1','DeltaM2'])
plot2d(mutArr1)
print(df)
df.to_excel('test1.xlsx', sheet_name='sheet1')

