import numpy as np
import scipy.optimize
from scipy.optimize import minimize
from sklearn.preprocessing import scale,StandardScaler,MinMaxScaler
from scipy.optimize import minimize_scalar
from chainer import FunctionSet, Variable, optimizers, serializers,iterators,training,Chain
import chainer.functions as F
import chainer.links as L
import pickle

AminoAcide='GAVILPSTCMDNEQKRHFYW'
AminoAcide={AminoAcide[i]:i for i in range(len(AminoAcide))}

with open('min_max_L_and_VJ_dict.p','rb') as file:
    minL, maxL, VJ_dict, AA=pickle.load(file)
print('VJ',len(VJ_dict))
print('AA',len(AA))


class scipyModel:
    def __init__(self,DATA,GEN,count):
        self.QL=np.random.uniform(0.5,1.5,int(maxL-minL+1))
        self.QVJ = np.random.uniform(0.5, 1.5, len(VJ_dict))
        self.QA = np.random.uniform(0.5, 1.5, len(AA)+1)
        self.QA[0]=1
        self.Z=1

        self.l_GEN_indexes = {}
        self.l_DATA_indexes = {}
        self.vj_GEN_indexes = {}
        self.vj_DATA_indexes = {}
        self.aa_GEN_indexes = {}
        self.aa_DATA_indexes = {}
        self.GEN = GEN
        self.DATA = DATA
        self.count=count
        self.S=np.sum(self.count)
        self.used_coef = np.unique(DATA[:, 2:])

        self.AA_back_code={val:key for key,val in AA.items()}

        print(1)
        for k in range(maxL-minL+1):
            self.l_DATA_indexes[k] = np.sum(self.count[np.where(DATA[:, 0] == k)[0]])
            self.l_GEN_indexes[k] = np.where(GEN[:, 0] == k)[0]
        print(2)
        for k in VJ_dict.values():
            self.vj_DATA_indexes[k] = np.sum(self.count[np.where(DATA[:, 1] == k)[0]])
            self.vj_GEN_indexes[k] = np.where(GEN[:, 1] == k)[0]
        print(3)
        for k in AA.values():
            column=self.AA_back_code[k][1]
            self.aa_DATA_indexes[k] = np.sum(self.count[np.where(DATA[:, 2+column] == k)[0]])
            self.aa_GEN_indexes[k] = np.where(GEN[:, column+2] == k)[0]

    def get_AA_coeficient(self):
        out = np.zeros((20, maxL - minL + 1, maxL))
        for i in range(1, len(self.QA)):
            if i in self.used_coef:
                out[AminoAcide[self.AA_back_code[i][2]], self.AA_back_code[i][0] - minL, self.AA_back_code[i][1]] = self.QA[i]
        return out

    def predict(self,DATA):
        o1 = self.QL[DATA[:, 0]]
        o2 = self.QVJ[DATA[:, 1]]
        shape = DATA[:, 2:].shape
        o3 = self.QA[list(DATA[:, 2:].reshape(-1))].reshape(shape)
        return np.prod(o3, axis=1) * o2 * o1/self.Z
    def evalf_Q(self,x,DATA):
        QL = x[:int(maxL - minL + 1)]
        QVJ = x[int(maxL - minL + 1):int(maxL - minL + 1) + len(VJ_dict)]
        QA = np.array([1] + list(x[int(maxL - minL + 1) + len(VJ_dict):]))
        o1 = QL[DATA[:, 0]]
        o2 = QVJ[DATA[:, 1]]
        shape = DATA[:, 2:].shape
        o3 = QA[list(DATA[:, 2:].reshape(-1))].reshape(shape)
        return np.prod(o3,axis=1)*o2*o1

    def likehood(self,x):
        QL = x[:int(maxL - minL + 1)]
        QVJ = x[int(maxL - minL + 1):int(maxL - minL + 1) + len(VJ_dict)]
        QA = np.array([1] + list(x[int(maxL - minL + 1) + len(VJ_dict):]))
        o1 = QL[list(self.DATA[:, 0])]
        o2 = QVJ[list(self.DATA[:, 1])]
        shape = self.DATA[:, 2:].shape
        o3 = QA[list(self.DATA[:, 2:].reshape(-1))].reshape(shape)

        Z=np.mean(self.evalf_Q(x,self.GEN))
        out=np.sum(self.count*(np.log(o1)+np.log(o2)+np.log(np.prod(o3,axis=1))-np.log(Z)))
        print(out)
        return -out

    def gradient(self,x):
        QL = x[:int(maxL - minL + 1)]
        QVJ = x[int(maxL - minL + 1):int(maxL - minL + 1) + len(VJ_dict)]
        QA = np.array([1] + list(x[int(maxL - minL + 1) + len(VJ_dict):]))
        out=[]
        Z = np.sum(self.evalf_Q(x, self.GEN))
        for k in range(len(QL)):
            out.append((self.l_DATA_indexes[k]-self.S*np.sum(self.evalf_Q(x,self.GEN[self.l_GEN_indexes[k]])) / Z)/QL[k])
        for k in range(len(QVJ)):
            out.append((self.vj_DATA_indexes[k] - self.S * np.sum(self.evalf_Q(x, self.GEN[self.vj_GEN_indexes[k]])) / Z)/QVJ[k])
        for k in range(1,len(QA)):
            out.append((self.aa_DATA_indexes[k] - self.S * np.sum(self.evalf_Q(x, self.GEN[self.aa_GEN_indexes[k]])) / Z)/QA[k])

        return -np.array(out)

    def x0(self):
        return np.array(list(self.QL)+list(self.QVJ)+list(self.QA[1:]))

class Model:
    def __init__(self):
        pass
    def predict(self,DATA):
        return self.M.predict(DATA)

    def get_AA_coeficient(self):
        return self.M.get_AA_coeficient()

    def fit(self,DATA,GEN,count=None):
        if count==None:
            count=np.ones((len(DATA)))
        self.M=scipyModel(DATA,GEN,count)
        bounds = [(0.001, None) for x in self.M.x0()]
        res = minimize(self.M.likehood, self.M.x0(), jac=self.M.gradient, method='L-BFGS-B', bounds=bounds,options={'maxiter':500})
        QL = res.x[:int(maxL - minL + 1)]
        QVJ = res.x[int(maxL - minL + 1):int(maxL - minL + 1) + len(VJ_dict)]
        QA = res.x[int(maxL - minL + 1) + len(VJ_dict):]
        self.M.QL=QL
        self.M.QVJ=QVJ
        self.M.QA=np.array([1]+list(QA))
        self.M.Z=np.mean(self.M.evalf_Q(res.x,GEN))







class OldModel():
    def __init__(self):
        self.model = Chain(
            QL=L.EmbedID(int(maxL-minL+1),1),
            QVJ=L.EmbedID(len(VJ_dict), 1),
            QA=L.EmbedID(len(AA)+1,1,ignore_label=0)

        )

        self.optimizer = optimizers.Adam()
        self.optimizer.setup(self.model)
        self.Z=1


    def forward(self,data):

        out1=self.model.QL(data[:,[0]])
        out2 = self.model.QVJ(data[:, [1]])
        out3=self.model.QA(data[:,2:])
        out=F.sum(F.concat((out1,out2,out3),axis=1)[:,:,0],axis=1)
        return out

    def get_AA_coeficient(self):
        QA=self.model.QA.W.data
        back_code={val:key for key,val in AA.items()}
        out=np.zeros((20,maxL-minL+1,maxL))
        for i in range(1,len(QA)):
            if i in self.used_coef:
                out[AminoAcide[back_code[i][2]],back_code[i][0]-minL,back_code[i][1]]=np.exp(QA[i])
        return out

    def predict(self,x):
        print('prediction')
        out=[]
        step=50
        for i in range(0,len(x),step):
            px = Variable(x[i:i+step])
            output = F.exp(self.forward(px))/self.Z
            out+=list(output.data)
        return np.array(out)



    def fit(self,DATA,GEN,counts,epochs):
        # Make a training function which takes in training data and targets
        # as an input.

        def train(data, gen, model,counts, batchsize=200000, n_epochs=1):
            lcurve=[]
            data_size = data.shape[0]
            gen_size = gen.shape[0]
            gen_var = Variable(gen)

            #self.model.to_gpu()
            for epoch in range(n_epochs):


                # randomly shuffle the indices of the training data
                shuffler = np.random.permutation(data_size)

                # loop over batches

                for i in range(0, data_size, batchsize):
                    print('epoch %d %d%%' % (epoch + 1,100*i/data_size))
                    data_var = Variable(data[shuffler[i: i + batchsize]])
                    counts_var = Variable(counts[shuffler[i: i + batchsize]])
                    S=np.sum(counts[shuffler[i: i + batchsize]])
                    output_d = self.forward(data_var)
                    output_g = F.exp(self.forward(gen_var))

                    model.zerograds()

                    loss = -(F.sum(output_d*counts_var)-S*F.log(F.sum(output_g)/gen_size))
                    lcurve.append(loss.data)
                    loss.backward()
                    self.optimizer.update()
            return np.array(lcurve)


        self.training=True
        self.lcurve=train(DATA, GEN, self.model, counts, n_epochs=epochs)

        gen_size = GEN.shape[0]
        gen_var = Variable(GEN)
        output_g = F.exp(self.forward(gen_var))
        self.Z = F.sum(output_g) / gen_size

        self.training=False
        self.used_coef=np.unique(DATA[:,2:])



