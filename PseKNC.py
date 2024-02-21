import math
import numpy as np

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


#فراوانی 2مرها رو محاسبه میکنه
def f2_kmer(SSeq):
    Nuc2=["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]
    num = 0
    all_freq2 = []
    for seq in SSeq:
        f2 = []
        for N in Nuc2:
            for i in range(len(seq) - 1):
                if seq[i:i+2]==N:
                    num=num+1
                else:
                    continue
            f2.append(num/199)
            num = 0
        all_freq2.append(f2)
    # print(all_freq2)
    return all_freq2


# بدست آوردن فرمول تتا
def correlation_factor_PseDNC(SSeq,SCValues,iter):
    list=[]
    for seq in SSeq:
        Theta = []
        for j in range(1, iter + 1):
            SUM = []
            for i in range(0, len(seq)-j-1):
                R1 = seq[i:i+2]
                R2 = seq[i+j:i+j+2]
                index1 = SCValues.index(R1)
                index2 = SCValues.index(R2)
                sum=0
                for z in range(1, 39):
                    PC1=float(SCValues[index1 + z])
                    PC2=float(SCValues[index2 + z])
                    sum+=math.pow((PC1-PC2),2)
                sum=sum/38
                SUM.append(sum)
            SUM=np.array(SUM)
            theta=np.sum(SUM)/(199-j)
            Theta.append(theta)
        list.append(Theta)
    return list

#بدست آوردن بردار du
def PseDNC(corfactor,freq,iter,w):
    Nuc2 = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
    list=[]
    Vetor = []
    for fr,cor in zip(freq,corfactor):
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


#فراوانی 3مرها رو محاسبه میکنه
def f3_kmer(SSeq):
    Nuc3=["AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
          "CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
          "GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
          "TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"]
    num = 0
    all_freq3 = []
    for seq in SSeq:
        f3 = []
        for N in Nuc3:
            for i in range(len(seq) - 2):
                if seq[i:i+3]==N:
                    num=num+1
                else:
                    continue
            f3.append(num/198)
            num = 0
        all_freq3.append(f3)
    return all_freq3

# بدست آوردن فرمول تتا
def correlation_factor_PseTNC(SSeq,SCValues,iter):
    list=[]
    for seq in SSeq:
        Theta = []
        for j in range(1, iter + 1):
            SUM = []
            for i in range(0, len(seq)-j-2):
                R1 = seq[i:i+3]
                R2 = seq[i+j:i+j+3]
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
def PseTNC(corfactor,freq,iter,w):
    Nuc3=["AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
          "CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
          "GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
          "TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"]
    Vetor = []
    for fr,cor in zip(freq,corfactor):
        vet = []
        fr = np.array(fr)
        Freq = np.sum(fr)
        cor = np.array(cor)
        Cor = np.sum(cor)
        for i in range(len(Nuc3)):
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


def Data_process_PseDNC(SSeq,valuesfilename,iter,w):
    Freq= f2_kmer(SSeq)
    valueslist = getSCValues(valuesfilename)
    corfactor= correlation_factor_PseDNC(SSeq, valueslist, iter)
    feature = PseDNC(corfactor, Freq, iter, w)

    Data= np.array(feature)

    return Data

def Data_process_PseTNC(SSeq,valuesfilename,iter,w):
    Freq= f3_kmer(SSeq)
    valueslist = getSCValues(valuesfilename)
    corfactor= correlation_factor_PseTNC(SSeq, valueslist, iter)
    feature = PseTNC(corfactor, Freq, iter, w)

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



best_parameter=best_parameter = {'k':2, 'iter': 3,'w': 0.2}
Sseq12,Label12=load_data(filename1,filename2)

print("hello")
if(best_parameter['k']==2):
    Data12 = Data_process_PseDNC(Sseq12, valuesfile_D, best_parameter['iter'], best_parameter['w'])
    np.savetxt('PseDNC', Data12)

elif(best_parameter['k']==3):
    Data12 = Data_process_PseTNC(Sseq12, valuesfile_T, best_parameter['iter'], best_parameter['w'])
    np.savetxt('PseTNC', Data12)

else:
    print("error2")