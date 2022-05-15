import numpy as np
import math
import random
import weblogo as wl

#Processing file
def getSeq(file):
    fh=open(file)
    myfile=fh.readlines()
    fh.close()
    seq=[]
    for i in range(1,len(myfile),2):
        seq.append(myfile[i].strip())
    return seq

def PosFreqMat(seq):
    """It will take a motif sequence and return the frequency of each nucleotide at position j ranging from 1 to length of sequence"""
    pfm=np.zeros((4,len(seq[0])))
    for i in range(len(seq)):
        for j in range(len(seq[0])):
                if seq[i][j]=="A":
                    pfm[0][j]=pfm[0][j]+1
                elif seq[i][j]=="U":
                    pfm[1][j]=pfm[1][j]+1
                elif seq[i][j]=="G":
                    pfm[2][j]=pfm[2][j]+1
                elif seq[i][j]=="C":
                    pfm[3][j]=pfm[3][j]+1          
    return pfm

def PosWtMat(pfmMatrix,seq):
    """It will take position frequency matrix and motif sequence and return the weight of each nucleotide at position j ranging from 1 to length of sequence"""
    pwm=np.zeros((4,len(seq[0])))
    #pwmlog=numpy.zeros((4,len(seq[0])))
    for i in range(4):
        for j in range(len(pfmMatrix[0])):
            pwm[i][j]=((pfmMatrix[i][j]+0.5)/(len(seq)+2)).round(4)
            #pwmlog[i][j]=math.log((((pfmMatrix[i][j]+0.5)/(len(seq)+2)).round(4))/(0.25))
    return pwm

def logoutscore(pwmMatrix,longseq,kmer):
    """ Calculating the weight of all possible window of random sequence and normalizing their value """
    lir=[]
    for i in range(len(longseq)-kmer+1):
        w=1.0
        subseq=longseq[i:i+kmer]
        for k in range(len(subseq)):
                if subseq[k]=="A":
                    w=w*pwmMatrix[0][k]
                elif subseq[k]=="U":
                    w=w*pwmMatrix[1][k]
                elif subseq[k]=="G":
                    w=w*pwmMatrix[2][k]
                elif subseq[k]=="C":
                    w=w*pwmMatrix[3][k]
        
        logw=math.log(w/((0.25)**len(subseq)))
        lir.append(logw)
    return lir

def getMotif(seq, kmer, sims, seed):
    #Randomly initializing the position motif in each sequence
    motifpos=[0]*len(seq)
    motifseq=[0]*len(seq)
    for i in range(len(seq)):
        motifpos[i]=random.choice(range(len(seq[i])-kmer+1))
        motifseq[i]=seq[i][motifpos[i]:motifpos[i]+kmer]

    #Iterating it for several times
    random.seed(seed)
    for lp in range(sims):
        #randomly select one sequence
        rs=random.choice(range(len(seq)))
        z=seq[rs]
        del motifseq[rs]
        del motifpos[rs]
        pfmMatrix=PosFreqMat(motifseq)
        pwmMatrix=PosWtMat(pfmMatrix,motifseq)
        y=logoutscore(pwmMatrix,z,kmer)

        windowS=[]
        for l in y:
            if l>0:
                windowS.append(l)
            else:
                windowS.append(0)

        mx=max(windowS)
        for i in range(len(windowS)):
            windowS[i]=windowS[i]/mx

        sum=0
        rp=random.uniform(0,1)
        for i in range(len(windowS)):
            sum=sum+windowS[i] 
            if sum>rp:
                cWin=i
                break;

        motifpos.insert(rs,cWin)
        motifseq.insert(rs,seq[rs][cWin:cWin+kmer])
    return motifpos,motifseq

def getEPS(motifpos, motifseq):
    #Writing fasta for motifseq as requiref by weblogo
    ofile = open("motifs.fa", "w")
    for i in range(len(motifseq)):
        ofile.write(">" + str(motifpos[i]) + "\n" +motifseq[i] + "\n")
    ofile.close()
    
    #getting eps from weblogo
    fin = open('motifs.fa')
    seqs = wl.read_seq_data(fin)
    logodata = wl.LogoData.from_seqs(seqs)
    logooptions = wl.LogoOptions()
    logooptions.title = "Sequence Logo"
    logoformat = wl.LogoFormat(logodata, logooptions)
    eps = wl.eps_formatter(logodata, logoformat)
    return eps

def runGS(file,kmer,sims,seed):
    seq=getSeq(file)
    mp,ms=getMotif(seq, kmer, sims, seed)
    eps=getEPS(mp,ms)

    with open(file[:-3]+"_"+str(seed)+"seed_"+str(kmer)+"mer_"+str(sims)+"sims.eps", "wb") as fo:
        fo.write(eps)
    pass