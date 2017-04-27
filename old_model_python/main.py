import numpy as np
import pickle
import matplotlib.pyplot as pl
from model import Model


AA='GAVILPSTCMDNEQKRHFYW'




DATA=np.load('X_data.npy')
GEN=np.load('X_gen.npy')
counts=np.array(np.load('counts.npy'),dtype=np.float32)

print(DATA.shape,GEN.shape,counts.shape)


M=Model()




M.fit(DATA,GEN,counts)
#exit()

with open('Y_model2.p','wb') as file:
    pickle.dump(M,file)


QA=M.get_AA_coeficient()
np.save('QA2.npy',QA)

"""pl.figure()
pl.plot(M.lcurve)

for i in range(5):
    pl.figure()
    pl.title(AA[i])
    pl.imshow(QA[i],interpolation='none')"""

#np.save('lcurve.npy',M.lcurve)
