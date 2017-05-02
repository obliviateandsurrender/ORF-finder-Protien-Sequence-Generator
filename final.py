import csv
import re
import textwrap
from operator import add
from functools import reduce

in_file1 = open("protein.fasta","r")
in_file2 = open("Introns.csv")
op1_file = open("1.txt","w")
op2_file = open("2.txt","w")
op3_file = open("3.txt","w")
op4_file = open("4.txt","w")
op5_file = open("5.txt","w")

pattern = re.compile(r'(ATG(?:[ATGC]{3})*?([ATGCZ*])*?(?:TAA|TAG|TGA))') #ORF SEUQENCE 
pattern1 = re.compile(r'([CG]{4,30})')

def flatten_List (matrix): #MAKES 2D LIST A 1D LIST
    vector = reduce(add, matrix, [])
    return vector

def ReverseComplement1(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','*':'*','Z':'Z'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def orfs(ira): #USE THE PATTER TO FIND ANSWER OF MY LIFE
    return list(pattern.findall(ira))

def translate_dna(sequence):

    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
    'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    proteinsequence = ''

    for n in range(0,len(sequence),3):
        if sequence[n:n+3] in codontable:
            proteinsequence += codontable[sequence[n:n+3]]
        else:   
            proteinsequence += ''
    return proteinsequence

def intron_remover(sequence):
    for i in sequence:
        l3 = []
        k = ''
        if (int(s2[0])-1 > 0):
            k = i[0:int(s2[0])-1]
        for j in range(0, len(s2),2):
            #print(j)
            for l in range(int(s2[j]),int(s2[j+1])+1):
                #print(int(s2[j]))
                #print(int(s2[j+1])+1)
                k += '*'
            if (j+2 < len(s2)):
                k += i[int(s2[j+1]):int(s2[j+2])-1]
            #print (len(k))
            if (j+2 >= len(s2)):
                break
        #    print(j)
        if (int(s2[-1]) < len(i)):
            k += i[int(s2[j-1]):-1]
            k += i[-1]
        #print (i)
        return k

def atgc(sequence):
    freq = [0,0,0,0]
    for i in sequence:
        if i == 'A':
            freq[0]+=1
        elif i == 'T':
            freq[1]+=1
        elif i == 'G':
            freq[2]+=1
        elif i == 'C':
            freq[3]+=1
        else:
            continue
    return freq

def tm(freq):
    return ( 2*(freq[0]+freq[1]) + 4*(freq[2]+freq[3]))
    

def frek(freq):
    ratio = (freq[2]+freq[3]) / (freq[0] + freq[1] + freq[2] + freq[3])
    if (ratio >= .40 and ratio <= .60):
        return 1
    else:
        return 0

def clust(sequence):
    if (pattern1.match(sequence)):
        #print(pattern1.findall(sequence))
        return 0
    else:
        return 1

def lrun(sequence):
    pattern2 = re.compile(r'([A]{5,30})')
    if (pattern2.match(sequence)):
        #print(pattern2.match(sequence))
        return 0
    pattern2 = re.compile(r'([T]{5,30})')
    if (pattern2.match(sequence)):
        #print(pattern2.match(sequence))
        return 0
    pattern2 = re.compile(r'([G]{5,30})')
    if (pattern2.match(sequence)):
        #print(pattern2.match(sequence))
        return 0
    pattern2 = re.compile(r'([C]{5,30})')
    if (pattern2.match(sequence)):
        #print(pattern2.match(sequence))
        return 0
    return 1

def wobble(sequence):
    if (sequence[-1] == 'A' or sequence[-1] == 'C'):
        return 1
    else:
        return 0

def end3(freq):
    ratio = (freq[1]+freq[0]) / (freq[0] + freq[1] + freq[2] + freq[3])
    if (ratio >= .50):
        return 1
    else:
        return 0

def hair(sequence1,sequence2):
    sequence2 = (ReverseComplement1(sequence2))[::-1]
    if (sequence1 == sequence2):
        return 0
    else:
        return 1

def primer(sequence,lmno):
    l = []
    flag = 0
    cunt = 0
    #print (sequence)
    if ((lmno+1) < 4): 
        slt = "> 5`-3` ORF" + str(counter%4+1)        
    else:
        slt = "> 3`-5` ORF" + str((counter+2)%4)
    op4_file.write(slt)
    op4_file.write("\n")
    for j in range(18,30):
        flag = 0
        for i in range(0, len(sequence)):
            sum = 0
            if (j+i > len(sequence)):
                break

            seq = (sequence)[i:j+i]
            user = atgc(seq)
            t = tm(user)
            if (t > 60 and t < 70):
                sum +=1
            sum += frek(user)
            sum += clust(seq)
            sum += lrun(seq)
            sum += wobble(seq)
            sum += end3(user)
            sum += hair(seq[0:9],seq[-10:-1])

            if (sum == 7):
                l.append(seq)
                cunt +=1

            if (cunt == lmno+2):
                flag = 1
                break
        
        if (flag):
            break
    flag = 0
    #print (l)
    
    k = l[lmno]
    #print (k)
    for j in range(18,30):
        flag = 0
        for i in range(0, len(sequence)):
            sum = 0
            seq = (ReverseComplement1(sequence))[i:j+i]
            if (j+i > len(sequence)):
                break
            user = atgc(seq)
            t = tm(user)
            if (t > 60 and t < 70):
                sum +=1
            sum += frek(user)
            sum += clust(seq)
            sum += lrun(seq)
            sum += wobble(seq)
            sum += end3(user)
            sum += hair(seq[0:9],seq[-10:-1])

            jpg = tm(atgc(k))
            if (sum ==7  and abs(jpg-t) < 5 and abs(jpg-t) > lmno%4):
                string = "Forward Primer - " + str(k) + "\t" + "Tm(in deg C) - " + str(jpg)
                #print (string)
                op4_file.write(string)
                op4_file.write("\n")
                string = "Reverse Primer - " + str(seq) + "\t" + "Tm(in deg C) - " + str(tm(user))
                #print (string)
                op4_file.write(string)
                op4_file.write("\n")
                flag = 1
                break
        if (flag):
            break
    return 


s1 = in_file1.read()
s2 = list(csv.reader(in_file2))
s2.pop(0)
s2 = flatten_List(s2)
 
l1 = s1.split(">")
l1.pop(0)

l2=['*']

info = []
sequence = []
isequence = []

#print(s2)

for i in l1:
    l3 = i.split("\n")
    info.append(l3[0])
    l3.pop(0)
    sequence.append(''.join(l3))


    #        flag = 1
    #        break
    #if (flag):
    #    break

isequence.append(intron_remover(sequence))
#print (isequence[0])
orf1 = []
iorf1 = []
#print (isequence)
#print (sequence[0])
#print (ReverseComplement1(sequence[0]))

for i in range(0,3):
    orf1.append(sequence[0][i:]) 
    iorf1.append(isequence[0][i:])
#print (sequence[0])

revsequence = ReverseComplement1(sequence[0]) #[::-1]
irevsequence = ReverseComplement1(isequence[0])
#print (revsequence)
#print (reverse(irevsequence))

for i in range(0,3):
    orf1.append(revsequence[i:])
    iorf1.append(irevsequence[i:])

orf2 = []
orf7 = []
counter = 0

for i in orf1:
    primer(str(i),counter)
    counter+=1

for i in iorf1:
    stseq = ''
    count = 0
    for j in i:
        count+=1
        stseq += j 
        if (count%3 == 0):
            stseq += 'Z'
    orf2.append(stseq)

pmer = (ReverseComplement1(orf2[0])[::-1])
pmer = pmer.replace("*","")
pmer = pmer.replace("Z","")
#print (textwrap.fill(pmer))
for i in orf1:
    stseq = ''
    count = 0
    for j in i:
        count+=1
        stseq += j 
        if (count%3 == 0):
            stseq += 'Z'
    orf7.append(stseq)

#a = orf1[0].replace("Z","")

#print (a)

counter = 0
for i in orf7:
    counter += 1
    if (counter < 4): 
        slt = "> 5`-3` ORF" + str(counter%4)
    else:
        slt = "> 3`-5` ORF" + str((counter+1)%4)
    op1_file.write(slt)
    op1_file.write("\n")
    op1_file.write(textwrap.fill(i.replace("Z",""),70))
    op1_file.write("\n")
    op1_file.write("\n")


counter = 0
length = 0
cunt  = 0
fsquirt = ''

for i in orf2:
    a = i.replace("Z"," ")
    a = a.replace("*","X")
    #print (a)
    #print (counter+1)
    for key1 in range(0,len(i),4):
        if (i[key1:key1+3] == 'ATG'):
            break;
    for key2 in range(len(i)-1,-1,-1):
        if (i[key2] != '*' and i[key2] != 'Z'):
            break;
    i = i[key1:key2+1]
    a = i.replace("Z","")
    a = a.replace("*","")
    #print (textwrap.fill(a,70))
    stseq = ''
    count = 0
    for mem in a:
        count+=1
        stseq += mem 
        if (count%3 == 0):
            stseq += 'Z'
    a = stseq
    a = orfs(a)
    counter += 1
    #print (counter)
    lona = 0
    for j in a:
        j = list(j)
        j.pop(-1)
        j[0] = j[0].replace("Z","")
        j[0] = j[0].replace("*","X")
        #print (len(j[0]))
        #print (textwrap.fill(j[0]))
        if (length < len(j[0])):
            cunt  = counter
            fsquirt = j[0]
            length = len(j[0])

if (cunt < 4): 
    slr = "> 5`-3` ORF" + str(cunt%4)
else:
    slr = "> 3`-5` ORF" + str((cunt+1)%4)
op2_file.write(slr)
op2_file.write("\n")
op2_file.write(textwrap.fill(orf7[cunt-1].replace("Z",""),70))


op3_file.write("> Protein Sequence GH")
op3_file.write("\n") 
op3_file.write(textwrap.fill(translate_dna(fsquirt),70))

op5_file.write("> E.coli cloning sequence")
op5_file.write("\n")
op5_file.write(textwrap.fill(fsquirt,70))
op5_file.write("\n")
op5_file.write("\n")
op5_file.write("> Pichia pastoris cloning sequence")
op5_file.write("\n")
op5_file.write(textwrap.fill(orf7[cunt-1].replace("Z",""),70))
op5_file.write("\n")
op5_file.write("\n")
op5_file.write("> HEK293 cloning sequence")
op5_file.write("\n")
op5_file.write(textwrap.fill(orf7[cunt-1].replace("Z",""),70))
op5_file.write("\n")
op5_file.write("\n")

in_file1.close()
in_file2.close()
