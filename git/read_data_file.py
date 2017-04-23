import numpy as np
from matplotlib import pyplot as pl

seq=[]
p=[]
y=[]
K=0

pgen=[]
with open('newdata/yzh.txt.ymir_probs_aa.txt') as file:
    for row in file:
        pgen.append(float(row[:-1]))
pgen=np.array(pgen)
pgen_out=[]

with open('newdata/yzh.txt') as file:
    print(len(next(file).split('\t')))
    for row in file:
        r=row.split('\t')
        K+=1
        if len(r)==16:
            if r[5].isalpha() and len(r[4])==len(r[5])*3 and pgen[K]>0:
                seq.append((r[5],r[6],r[7]))
                if '"' in r[7]:
                    print(';')
                p.append(float(r[0]))
                pgen_out.append(pgen[K])


pgen_out=np.array(pgen_out)

print(len(pgen_out),len(seq))

seq=np.array(seq)
print(seq.shape)
np.save('Y_data.npy',seq)
np.save('p_exp.npy',p)
np.save('p_gen.npy',pgen_out)




