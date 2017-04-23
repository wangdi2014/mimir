import numpy as np
from matplotlib import pyplot as pl

seq=[]
p=[]
y=[]
K=0



with open('newdata/yzh.aa.txt') as file:
    print(len(next(file).split('\t')))
    for row in file:
        r=row.split('\t')
        K+=1
        if len(r)==3:
            if r[0].isalpha():
                seq.append((r[0],r[1],r[2].replace('\n','')))
seq=np.array(seq)
print(seq.shape)
np.save('Y_gen.npy',seq)




