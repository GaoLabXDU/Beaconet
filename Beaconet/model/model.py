# it is an implementation of Beaconet
import numpy as np
from tqdm import trange
import torch as t
from torch import nn,autograd
from torch.optim import Adam
from functools import reduce
import pandas as pd
import warnings

"""
Beaconet is a reference-free method for integrating scRNA-seq datasets. Its major properties are reference-free and
working in the original molecular feature space. The reference-free property is mainly for overcoming the reliance of the
pre-selected reference for integration (different reference may lead to different integration performance due to the
cell diversity and the data quality of batches). The original feature space preserved working fashion is mainly designed
for enhancing more potential downstream analysis which not only need identify the structure of data (such as cluster)
but also the molecular features without batch effect. 

For conveniently calling, we provide a top level encapsulation of the implementation of Beaconet.
***correction***
"""


def correction(dfs,device=None,n_critic=10,Lambda=10,d_model=256,minibatch_size=1024):
    """
        This is the API to calling Beaconet for integrating scRNA-seq datasets.
        see more details in our paper.

        Return:
        ----------
        res: pandas.DataFrame, the integrated dataset of the given batches of scRNA-seq data.

        Parameters
        ----------
        dfs  : List of pandas.DataFrame, each element of the list is a log-scaled transcriptomic dataset,
            the shape of each matrix is (n_cells, n_genes).
            The number of columns (gene) must be equal among matrices.
            The number of rows (cell) may be different among matrices.

        device: None or device. indicate the device for computing. default None, the default behaviour is to use
            gpu 'cuda:0' if available, to use cpu otherwise.
            The user also can set device by passing the arguement using the return value of torch.device().

        n_critic: int, default 10. In the adversarial training process of Beaconet, the Corrector and Discriminator are
            trained alternately. the parameter *n_critic* controls the balance of the two components. Specificly, in each
            epoch, the Discriminator is trained in *n_critic* step while Corrector is trained in one step.

        Lambda: float, default 10. Lambda is the hyper-parameter for Lipthiz-1 condtion penalty term in loss of Beaconet.
            This term is originally presented in this paper:
            Gulrajani, I., Ahmed, F., Arjovsky, M., Dumoulin, V., & Courville, A. C. (2017). Improved training of wasserstein gans. Advances in neural information processing systems, 30.

        d_model: int, default 256. The number of dimensions of the latent space in Beaconet. For convenience and
            illustrating the effectiveness of Beaconet, we simply set the same number of dimensions in all fully connected
            layers.

        minibatch_size: int, default 1024. The size of mini-batch for stochastic gradient descent.

        =========
        The usage demo are included in the Beaconet/test/demo.py
        =========
    """
    if(device is None):
        device = t.device("cuda:0" if (t.cuda.is_available()) else "cpu")

    if(not t.cuda.is_available()):
        text="""
The cuda or gpu is not available. The calculation may be slow. We strongly suggest to apply Beacoent on GPU, since 
GPU and cuda can improve the efficiency and effectiveness significantly.
        """
        warnings.warn(text)

    ###################
    #train model
    model=Model(dfs,device,n_critic=n_critic,LAMBDA=Lambda,d_model=d_model,minibatch_size=minibatch_size)
    print("training...")
    model.train()
    print("calculating the results...")
    res = model.get_result(dataframe=True)
    print("finish!")
    return res


class Model:
    """
        This is the model class of Beaconet, it contains a Corrector and a Discriminator. The role of Corrector is to
        correct log-scaled scRNA-seq data. It accepts the vector of log-scaled gene-expression features and batch index as
        input and output the corrected scRNA-seq data. The role of Discriminator is to identify the batch difference among
        multiple datasets based on Wasserstein distance, and guide the optimization of Corrector.
        see more details in our paper.

        For the bioinformatic researcher who is not familiar to machine learning technologies, we suggest to use the
        top-level function "correction", which packaged the model, training strategy, objective function, and default
        value of several hyperparameters. "correction" function is more convenient for using.

    """
    def __init__(self,dfs,device,d_model=256,LAMBDA=10,n_critic=10,minibatch_size=1024):
        """
        Parameters
        ----------
        dfs  : List of pandas.DataFrame, each element of the list is a log-scaled transcriptomic dataset,
            the shape of each matrix is (n_cells, n_genes).
            The number of columns (gene) must be equal among matrices.
            The number of rows (cell) may be different among matrices.

        device: None or device. indicate the device for computing. default None, the default behaviour is to use
            gpu 'cuda:0' if available, to use cpu otherwise.
            The user also can set device by passing the arguement using the return value of torch.device().

        n_critic: int, default 10. In the adversarial training process of Beaconet, the Corrector and Discriminator are
            trained alternately. the parameter *n_critic* controls the balance of the two components. Specificly, in each
            epoch, the Discriminator is trained in *n_critic* step while Corrector is trained in one step.

        Lambda: float, default 10. Lambda is the hyper-parameter for Lipthiz-1 condtion penalty term in loss of Beaconet.
            This term is originally presented in this paper:
            Gulrajani, I., Ahmed, F., Arjovsky, M., Dumoulin, V., & Courville, A. C. (2017). Improved training of wasserstein gans. Advances in neural information processing systems, 30.

        d_model: int, default 256. The number of dimensions of the latent space in Beaconet. For convenience and
            illustrating the effectiveness of Beaconet, we simply set the same number of dimensions in all fully connected
            layers.

        minibatch_size: int, default 1024. The size of mini-batch for stochastic gradient descent.

        """
        n_batches = len(dfs)
        n_features=dfs[0].shape[1]
        assert all([d.shape[1]==n_features for d in dfs])
        self.total_cells = sum([e.shape[0] for e in dfs])
        self.C = Corrector(n_features=n_features, d_model=d_model, n_batches=n_batches).to(device=device)
        self.D = Discriminator(n_features=n_features, d_model=d_model, n_batches=n_batches).to(device=device)
        self.LAMBDA = LAMBDA
        self.n_critic=n_critic
        self.minibatch_size = minibatch_size
        self.data = Dataset(dfs, device, self.minibatch_size, norm=None)
        self.optC = Adam(self.C.parameters(), lr=2e-4, betas=(0, 0.999))
        self.optD = Adam(self.D.parameters(), lr=2e-4, betas=(0, 0.999))
        self.device=device
        self.iter=self.estimate_iter()

    def train(self,iter=None,iter_init=1000):
        if(iter is None):
            iter=self.iter
        else:
            self.iter=iter
        dataloader = self.data.inf_dataloader()
        print("initializing the W distance regression model")
        for i in trange(iter_init):
            minibatch = next(dataloader)
            self.updateD(minibatch=minibatch)

        print("training model")
        for i in trange(iter):
            # train D
            for j in range(self.n_critic):
                minibatch = next(dataloader)
                self.updateD(minibatch=minibatch)
            # train C
            self.updateC(minibatch=minibatch)


    def updateC(self,minibatch):
        D=self.D
        C=self.C
        optC=self.optC

        n_batches = len(minibatch)
        C.zero_grad()
        D.zero_grad()
        loss = []
        for i, (y, x) in enumerate(minibatch):
            x1, c = C(x, y)
            disc = D(x1).mean(dim=0)
            loss.append(disc)

        loss = t.stack(loss, dim=0)
        coeff = t.zeros_like(loss)
        coeff[:] = 1 / (n_batches - 1)
        coeff.fill_diagonal_(0)

        loss = (loss * coeff).sum(dim=1).mean()
        loss.backward()
        optC.step()


    def updateD(self,minibatch):
        # Sampling
        D=self.D
        C=self.C
        optD=self.optD

        D.zero_grad()

        n_batches = len(minibatch)
        disc = []
        x_corrected = []
        for i, (y, x) in enumerate(minibatch):
            with t.no_grad():
                x1, _ = C(x, y)
            x_corrected.append(x1)
            disc.append(D(x1).mean(dim=0))

        disc = t.stack(disc, dim=0)
        coeff = t.zeros_like(disc)
        coeff[:] = -1 / (n_batches - 1)
        coeff.fill_diagonal_(1)
        # disc is the matrix consist of E[Dj(C(xi))]. the coefficients of diag values was 1, others was -1/(n-1).

        loss = (disc * coeff).sum()
        penalty = calculate_gradient_penalty(D, x_corrected, device=self.device)
        loss += self.LAMBDA * penalty
        loss.backward()
        optD.step()

    def estimate_iter(self,max_iter=2500,a=20,b=2.8):
        return min(int(a*np.log10(self.total_cells)**b),max_iter)

    @t.no_grad()
    def get_result(self,dataframe=False):
        """

        :param dfs:
        :return:
        """
        res = []
        c_res = []

        for i in range(self.data.n_batches):
            c = []
            one_batch_res = []
            for y, x in self.data.dataloader_for_eval(data_index=i):
                x1, cc = self.C(x, y)
                one_batch_res.append(x1)
                c.append(cc)
            one_batch_res = t.cat(one_batch_res).cpu()
            c_res.append(t.cat(c).cpu())
            res.append(one_batch_res)

        if(dataframe):
            res = t.cat(res).numpy()
            return pd.DataFrame(res,index=reduce(lambda x,y:list(x)+list(y),self.data.cell_name))

        return res, c_res

    def save(self,filename=None):
        if(filename is not None):
            t.save(
                {
                    "C": self.C.state_dict(),
                    "D": self.D.state_dict(),
                    "optC": self.optC.state_dict(),
                    "optD": self.optD.state_dict()
                }, filename)

    def load(self,filename):
        state = t.load(filename)
        self.C.load_state_dict(state["C"])
        self.D.load_state_dict(state["D"])
        self.optC.load_state_dict(state["optC"])
        self.optD.load_state_dict(state["optD"])


class Discriminator(nn.Module):
    def __init__(self,n_features,d_model,n_batches):
        """
            Parameters
            ----------
            n_features  : int, the number of features, e.g., genes.

            d_model: int, default 256. The number of dimensions of the latent space in Beaconet. For convenience and
            illustrating the effectiveness of Beaconet, we simply set the same number of dimensions in all fully connected
            layers.

            n_batches: int, the number of batches, e.g., sequencing technologies.

        """
        super(Discriminator,self).__init__()
        self.n_features=n_features
        self.d_model=d_model
        self.n_batches=n_batches

        self.model=nn.Sequential(
            nn.Linear(n_features,d_model),
            nn.LeakyReLU(0.1),
            nn.Linear(d_model, d_model),
            nn.LeakyReLU(0.1),
            nn.Linear(d_model,n_batches),
        )

    def forward(self,x):
        return self.model(x)

class BatchSpecificNorm(nn.Module):
    def __init__(self,n_batches,n_features,eps=1e-8):
        super(BatchSpecificNorm, self).__init__()
        self.batch_c=nn.Embedding(n_batches,n_features)
        self.a = nn.Embedding(n_batches,1)

        #initialization of linear bias and scale
        self.batch_c.weight.data.zero_()
        #self.a.weight.data.fill_(eps)
        self.a.weight.data.fill_(1)

    def forward(self,x,y):
        return x * self.a(y) + self.batch_c(y)

class Corrector(nn.Module):
    def __init__(self,n_features,n_batches,d_model=512,eps=1e-8):
        """
        Parameters
        ----------
        n_features  : int, the number of features, e.g., genes.

        n_batches: int, the number of batches, e.g., sequencing technologies.

        d_model: int, default 256. The number of dimensions of the latent space in Beaconet. For convenience and
            illustrating the effectiveness of Beaconet, we simply set the same number of dimensions in all fully connected
            layers.

        """
        super(Corrector,self).__init__()

        self.fc=nn.Sequential(
            nn.Linear(n_features,d_model),
            nn.LeakyReLU(0.1),
            nn.Linear(d_model, d_model),
            nn.LeakyReLU(0.1),
            nn.Linear(d_model, n_features),
            nn.Tanh(),
        )
        self.bsn=BatchSpecificNorm(n_batches=n_batches,n_features=n_features,eps=eps)
        self.relu=nn.ReLU()

    def forward(self,x,y):
        c=self.bsn(self.fc(x),y)
        x_new=self.relu(x+c)# enforce the result expression of cells should not be non-negative values
        return x_new, x_new-x


def sampling_for_penalty(x_list):

    ms = np.array([x.shape[0] for x in x_list])
    ms=ms.max()-ms # how many extra sample was needed
    #n=len(x_list)
    new_x_list=[]
    for sub_m,x in zip(ms,x_list):
        x=x.detach()#
        x=t.cat([x,x[t.randint(high=x.shape[0], size=(sub_m,))]])
        new_x_list.append(x)

    return new_x_list

def construct_inter(x_list,index,device):
    """
    construct inter samples from x_list[index] and others.
    :param x_list:
    :param index:
    :param device:
    :return:
    """
    inter=[]
    cur=x_list[index]
    n=cur.shape[0]
    for i,x in enumerate(x_list):
        if(i!=index):
            eps = t.rand((n, 1), device=device)
            sub_inter= cur + eps * (x-cur)
            inter.append(sub_inter)

    inter=t.cat(inter)
    inter.requires_grad_()
    return inter

def calculate_gradient_penalty(net:nn.Module,x_corrected,device):
    #
    n_batches=len(x_corrected)
    x_list=sampling_for_penalty(x_corrected)


    penalty=[]
    for i in range(n_batches):
        inter = construct_inter(x_list, index=i, device=device)
        out = net(inter)
        out=out[:,i]
        gradient=autograd.grad(outputs=out,grad_outputs=t.ones_like(out),inputs=inter,retain_graph=True,create_graph=True, only_inputs=True)[0]
        slopes=t.sqrt(t.sum(gradient**2,dim=1))
        penalty.append( ((slopes-1)**2).mean() )

    return sum(penalty)#/len(penalty)


def numpy2tensor(df,device=None):
    return t.tensor(df, dtype=t.float32,device=device)

def datafram2tensor(df,device=None,norm=None):
    res=t.tensor(df.values, dtype=t.float32, device=device)
    if(norm is not None):
        #print(t.sum(res,dim=norm,keepdim=True).mean())
        res=res/t.sum(res,dim=norm,keepdim=True)
    return res

def random(n_samples=32,n_dim=64):
    """
    random variable generator. the output follows U[0,1).

    :param n_samples: int
    :param n_dim: int
    :return: torch.Tensor(n_samples,n_dim)
    """
    return t.randn(n_samples,n_dim)

class Dataset:
    def __init__(self,dfs,device,minibatch_size,norm):
        # the meaning of "batch" in batch_size and sub_batch_size was mini-batch
        # the "batch" in n_batches was the "batch" in the batch effect of single cell data
        self.cell_name = [e.index for e in dfs]
        self.device=device
        self.y_list = [t.tensor([i for _ in range(e.shape[0])], device=device) for i, e in enumerate(dfs)]
        self.data_list =[ datafram2tensor(e, device=device,norm=norm) for i, e in enumerate(dfs) ]
        self.minibatch_size=minibatch_size

        self.n_batches=len(self.y_list)

    def dataloader_for_eval(self,data_index):
        minibatch_size=self.minibatch_size
        data=self.data_list[data_index]
        y=self.y_list[data_index]
        for i in range(0, data.shape[0], minibatch_size):
            yield y[i:i + minibatch_size],data[i:i + minibatch_size]

    def inf_dataloader(self):
        """
        a dataloader that generates infinite data from Dataset.

        :return:
        """
        n=self.minibatch_size
        while(True):
            res=[(y[t.randint(high=data.shape[0],size=(n,))],
                  data[t.randint(high=data.shape[0], size=(n,))])
                for y,data in zip(self.y_list,self.data_list)]

            yield res