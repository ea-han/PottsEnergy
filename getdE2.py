#!/usr/bin/env python2
from scipy import *
import numpy as np
import os,sys,re
from io import StringIO


def alphabet_compare1(alphabet1):
    return "ABCD".index(alphabet1)*4

def alphabet_compare2(alphabet2):
    return "ABCD".index(alphabet2)

def other_alphabet(alphabet):
    return list("ABCD".replace(alphabet, ''))


def reduced_alphabet(reduxarray,pos,mut):
    pos_=int(pos)-1
    s=reduxarray[pos_]
    print(s)
    for k,txt in enumerate(s):
        if(txt.find(mut)!=-1):
            reduced_mut=" ABCD"[k]
    print(reduced_mut)
    return reduced_mut

seq=np.loadtxt(sys.argv[3],dtype='str')
redux=np.loadtxt(sys.argv[4],dtype='str')
J=np.load(sys.argv[2])
dummy,orig,p,mut,dummy=re.split(r'([A-Z]+)([0-9]+)([A-Z])',sys.argv[1])
pos=int(p)
print(pos, orig, mut)
reduced_orig=reduced_alphabet(redux,pos,orig)
reduced_mut=reduced_alphabet(redux,pos,mut)
print(sys.argv[1].lower())


Q_squared=J.shape[1]
Q=int(np.sqrt(Q_squared)+0.2)
L=int((1+np.sqrt(1+8*J.shape[0]))/2+0.2)
print(L,Q)
map_index=np.arange(0,int(J.shape[0]))
J_mapped=np.zeros([L,L])
J_mapped[np.triu_indices(L,1)]=map_index
J_mapped=J_mapped+J_mapped.T
print(J_mapped.shape)


conseq=np.loadtxt('in.consensus.reduce4.seq',dtype='str',ndmin=1)
con=list(conseq[0])
con2=[]
for c in con:
    con2.append(alphabet_compare2(c))


my_list1=[]
for i in range(0,len(seq)):
    posi=pos-1
    J_get=J_mapped[posi,:]
    gamma=con[posi]
    if(gamma!=reduced_orig):
        print("Consensus residue does not match!")
        break
    alpha=reduced_mut
    if(alpha==gamma):
            print("Mutation and consensus same in reduced alphabet! Stopping ...")
            break

    dE_tot=0
    for j,co in enumerate(J_get):
        J_ij=J[int(co)]        #or J_ij=J[int(J_get[j])]

        if(posi<j):
            E1=J_ij[(con2[posi]*Q)+alphabet_compare2(seq[i][j])]
            E2=J_ij[alphabet_compare1(alpha)+alphabet_compare2(seq[i][j])]
            dE=E1-E2
            #dE=E2-E1
        
        elif(posi==j):
            continue

        elif(posi>j):
            E1=J_ij[alphabet_compare1(seq[i][j])+con2[posi]]
            E2=J_ij[alphabet_compare1(seq[i][j])+alphabet_compare2(alpha)]
            dE=E1-E2
            #dE=E2-E1
        
        dE_tot=dE_tot+dE
    my_list1.append(dE_tot)

np.savetxt('dEs/de.exper.'+sys.argv[1].lower(),my_list1)
