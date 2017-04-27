import pickle
import numpy as np
X=np.load('Y_data.npy')
X_gen=np.load('Y_gen.npy')
p_exp=np.load('p_exp.npy')
#p_gen=np.load('p_gen.npy')


AA='GAVILPSTCMDNEQKRHFYW'


uniqL={}
for x,p in zip(X,p_exp):
    if len(x[0]) not in uniqL:
        uniqL[len(x[0])]=0
    uniqL[len(x[0])]+=int(p)
L=[]
for l,c in uniqL.items():
    if c>1000:
        L.append(l)
minL,maxL=np.min(L),np.max(L)

print(minL,maxL)
uV=[]
for u in X[:,1]:
    if ',' in u:
        uV.append(str(u).split(',')[0])
    else:
        uV.append(u)

uJ=[]
for u in X[:,2]:
    if ',' in u:
        uJ.append(str(u).split(',')[0])
    else:
        uJ.append(u)
uV=list(uV)+list(X_gen[:,1])
uJ=list(uJ)+list(X_gen[:,2])

VJ_dict={}


for v,j in zip(uV,uJ):
    if (v,j) not in VJ_dict:
        VJ_dict[(v,j)]=len(VJ_dict)


AA_dict={}
for l in range(minL,maxL+1):
    for i in range(l):
        for a in AA:
            AA_dict[(l,i,a)]=len(AA_dict)+1

with open('min_max_L_and_VJ_dict.p','wb') as file:
    pickle.dump((minL,maxL,VJ_dict,AA_dict),file)


fullSize=1+1+maxL
print(len(uV)*len(uJ))



print(len(X),fullSize)


p_exp=p_exp[np.where(np.array([len(x[0]) for x in X])>=minL)[0]]
#p_gen=p_gen[np.where(np.array([len(x[0]) for x in X])>=minL)[0]]
X=X[np.where(np.array([len(x[0]) for x in X])>=minL)[0]]
p_exp=p_exp[np.where(np.array([len(x[0]) for x in X])<=maxL)[0]]
#p_gen=p_gen[np.where(np.array([len(x[0]) for x in X])<=maxL)[0]]
X=X[np.where(np.array([len(x[0]) for x in X])<=maxL)[0]]

out=np.zeros((len(X),fullSize),dtype=np.int32)
print('start data')
for i in range(len(X)):
    out[i,0]=len(X[i,0])-minL
    v=X[i,1]
    if ',' in v:
        v=str(v).split(',')[0]
    j = X[i, 2]
    if ',' in j:
        j = str(j).split(',')[0]
    out[i,1]=VJ_dict[(v,j)]
    for j in range(len(X[i,0])):
        out[i,2+j]=AA_dict[len(X[i,0]),j,X[i,0][j]]
np.save('X_data.npy',out)
np.save('counts.npy',p_exp)
#np.save('data_gen_p.npy',p_gen)

print('start gen')


X=X_gen[np.where(np.array([len(x[0]) for x in X_gen])>=minL)[0]]
X=X[np.where(np.array([len(x[0]) for x in X])<=maxL)[0]]

out=np.zeros((len(X_gen),fullSize),dtype=np.int32)
for i in range(len(X)):
    out[i,0]=len(X[i,0])-minL
    v = X[i, 1]
    if ',' in v:
        v = str(v).split(',')[0]
    j = X[i, 2]
    if ',' in j:
        j = str(j).split(',')[0]
    out[i, 1] = VJ_dict[(v, j)]
    for j in range(len(X[i, 0])):
        out[i, 2 + j] = AA_dict[len(X[i, 0]), j, X[i, 0][j]]
np.save('X_gen.npy',out)
