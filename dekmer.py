import math
import numpy as np
from collections import Counter
from itertools import product
# using collections.OrderedDict.fromkeys()
from collections import OrderedDict


#گرفتن توالی ها
def getDNA(filename):
    fr = open(filename, 'r')
    SSeq=[]
    for line in fr.readlines():
        line = line.strip('\n')
        if(line[0]!='>'):
            str=line.upper()
            SSeq.append(str)
    return SSeq

#گرفتن خواص ساختار فیزیکوشیمیایی
def getSCValues(valuesfile):
    valueslist=[]
    fr = open(valuesfile, 'r')
    for lines in fr.readlines():
        lines = lines.strip('\n')
        lines=lines.split(",")
        for line in lines:
            valueslist.append(line)
    return valueslist

def split_into_overlapping_strings(input_string, k):
    overlapping_strings = []
    for i in range(len(input_string) - k + 1):
        substring = input_string[i:i+k]
        overlapping_strings.append(substring)
    return overlapping_strings

def f_kmer(SSeq,kmer,pos1):
     k = kmer
     s = ' '.join(SSeq)
     Nuc2 = ["".join(p) for p in product("ATGC", repeat=k)]
     Nuc3 = []
     for string in Nuc2:
        string = string[:pos1-1] + string[k - 1:k]
        Nuc3 += [string]
     res = []
     res = list(OrderedDict.fromkeys(Nuc3))
     Nuc4=split_into_overlapping_strings(s, k)
     Nuc5 = []
     for string2 in Nuc4:
       string2 = string2[:pos1-1] + string2[k - 1:k]
       Nuc5 += [string2]
     num = 0
     all_freq2 = []
     for N in Nuc5:
       f2 = []
       for j in res:
           for b in Nuc5:
             if b==j:
                num=num+1
             else:
                continue
           f2.append(num / 64)
           num = 0
       all_freq2.append(f2)
     return all_freq2


# بدست آوردن فرمول تتا
def correlation_factor_PseKNC(SSeq,SCValues,iter,kmer):
    list=[]
    for seq in SSeq:
        Theta = []
        for j in range(1, iter + 1):
            SUM = []
            for i in range(0, len(seq)-j-2):
                R1 = seq[i:i+kmer]
                R2 = seq[i+j:i+j+kmer]
                index1 = SCValues.index(R1)
                index2 = SCValues.index(R2)
                sum=0
                for z in range(1, 13):
                    PC1=float(SCValues[index1 + z])
                    PC2=float(SCValues[index2 + z])
                    sum+=math.pow((PC1-PC2),2)
                sum=sum/12
                SUM.append(sum)
            SUM=np.array(SUM)
            theta=np.sum(SUM)/(198-j)
            # print("ith theta",i,theta)
            Theta.append(theta)
        # print("jth Theta",j,Theta)
        list.append(Theta)
    # print("list.shape",np.array(list).shape)
    # print(list)
    return list
#بدست آوردن بردار du
def PseKNC(corfactor,freq,iter,w,k):

    Nuc2 = ["".join(p) for p in product("ATGC", repeat=k-1)]
    list=[]
    Vetor = []
    for fr, cor in zip(freq, corfactor):
        vet = []
        fr = np.array(fr)
        Freq = np.sum(fr)
        cor = np.array(cor)
        Cor = np.sum(cor)
        for i in range(len(Nuc2)):
            A = fr[i]
            B = float(Freq) + float(w * float(Cor))
            Du1 = A / B
            vet.append(Du1)
        for j in range(iter):
            A = w * float(cor[j])
            B = float(Freq) + float(w * float(Cor))
            Du2 = A / B
            vet.append(Du2)
        Vetor.append(vet)
    return Vetor


def Data_process_PseKNC(SSeq, valuesfilename, kmer, iter, w, pos1):
    Freq = f_kmer(SSeq,kmer,pos1)
    valueslist = getSCValues(valuesfilename)
    corfactor = correlation_factor_PseKNC(SSeq, valueslist, iter,kmer)
    feature = PseKNC(corfactor, Freq, iter, w,kmer)

    data_feature = np.array(feature)
    return data_feature



def label_Pos(feature):
    label_Pos=[]
    for i in range(len(feature)):
        label_Pos.append(1)
    return label_Pos

def label_Neg(feature):
    label_Neg=[]
    for i in range(len(feature)):
        label_Neg.append(0)
    return label_Neg


def load_data(filename1, filename2):
    Dnaseq1 = getDNA(filename1)
    Dnaseq2 = getDNA(filename2)

    data1_np = np.array(Dnaseq1)
    data2_np = np.array(Dnaseq2)
    data_12 = np.hstack([data1_np, data2_np])

    y1 = label_Pos(Dnaseq1)
    y2 = label_Neg(Dnaseq2)
    y1_np = np.array(y1)
    y2_np = np.array(y2)
    y12 = np.hstack([y1_np, y2_np])
    return data_12,y12







filename1="C:/Users/space/PycharmProjects/articles code/B_Enhancer.txt"
filename2="C:/Users/space/PycharmProjects/articles code/B_NonEnhancer.txt"
valuesfile_D = "C:/Users/space/PycharmProjects/articles code/dinucleotides_revise.csv"
valuesfile_T = "C:/Users/space/PycharmProjects/articles code//trinucleotides_revise.csv"



best_parameter=best_parameter = {'k':4, 'iter': 3,'w': 0.2,'pos':3}
Sseq12,Label12=load_data(filename1,filename2)

print("hello")
Data12 = Data_process_PseKNC(Sseq12, valuesfile_T, best_parameter['k'],best_parameter['iter'], best_parameter['w'],best_parameter['pos'])
np.savetxt('PseKNC', Data12)
