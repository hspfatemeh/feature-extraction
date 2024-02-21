import numpy as np

#برای محاسبه فراوانی ها دو قست شده بخاطر اینکه در محاسبه ی فراوانی مثبت ها نسبت به منفی ها تعداد توالی ها یکی بیشتر می باشد
def positiveside(Sseq12,Label12,i):
    Sseq122 = round(len(Sseq12)/2)
    Label122 = round(len(Label12)/2)
    occure = []
    i=i
    for j in range(0, 6):
        count = 0
        for k in range(0, Sseq122):
            if Sseq12[i][j] == Sseq12[k][j] and i != k:
                count = +1
            else:
                continue
        occure.append(count / Sseq122)
    return occure
def negativeside(Sseq12, Label12,i):
    Sseq122 = len(Sseq12)
    Label122 = len(Label12)
    occure = []
    i=i
    h=round(Sseq122/2)
    for j in range(0, 6):
            count = 0
            for k in range(h, Sseq122):
                if Sseq12[i][j] == Sseq12[k][j] and i != k:
                    count = +1
                else:
                    continue
            occure.append(count / ((Sseq122-h)+1))

    return occure

#برای محاسبه فراوانی منفی ها نسبت به مثبت ها باید یکی به تقسیم اضافه شود
def positiveside2(Sseq12,Label12,i):
    Sseq122 = round(len(Sseq12)/2)
    Label122 = round(len(Label12)/2)
    occure = []
    i=i
    for j in range(0, 6):
        count = 0
        for k in range(0, Sseq122):
            if Sseq12[i][j] == Sseq12[k][j] and i != k:
                count = +1
            else:
                continue
        occure.append(count / Sseq122+1)
    return occure


def negativeside2(Sseq12, Label12,i):
    Sseq122 = len(Sseq12)
    Label122 = len(Label12)
    occure = []
    i=i
    h = round(Sseq122 / 2)
    for j in range(0, 6):
            count = 0
            for k in range(h, Sseq122):
                if Sseq12[i][j] == Sseq12[k][j] and i != k:
                    count = +1
                else:
                    continue
            occure.append(count / 4)

    return occure



def Bi_profile(Sseq12,Label12):
    Sseq12=Sseq12
    Label12=Label12
    occure=[]
    occure2 = []
    for i in range(0,len(Sseq12)):
        if i<round(len(Sseq12)/2):
           f=positiveside(Sseq12, Label12,i)
           g=negativeside(Sseq12, Label12,i)
           occure.append(f+g)
           occure.append(Label12[i])
        else:
           f = positiveside2(Sseq12, Label12, i)
           g = negativeside2(Sseq12, Label12, i)
           occure.append(f + g)
           occure.append(Label12[i])

    return occure


def getDNA(filename):
    fr = open(filename, 'r')
    SSeq = []
    for line in fr.readlines():
         line = line.strip('\n')
         if (line[0] != '>'):
             str = line.upper()
             SSeq.append(str)
    return SSeq


def label_Pos(feature):
        label_Pos = []
        for i in range(len(feature)):
            label_Pos.append(1)
        return label_Pos


def label_Neg(feature):
        label_Neg = []
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
        return data_12, y12

filename1="C:/Users/space/PycharmProjects/articles code/B_Enhancer - Copy.txt"
filename2="C:/Users/space/PycharmProjects/articles code/B_NonEnhancer - Copy.txt"

Sseq12,Label12=load_data(filename1,filename2)
f=Bi_profile(Sseq12,Label12)
print(f)