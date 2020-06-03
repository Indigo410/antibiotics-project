import numpy as np
from scipy.stats import norm,pearsonr

class Simulate:
    def __init__(self,n,b1,b2,b3,vars=None):
        """
        @param n: number of samples
        @param b1: all 0 for direct mediation, some otu directly influenced target otu otherwise
        """
        # input variables
        self.n=n
        self.b1=b1
        self.b2=b2
        self.b3=b3

        # inferred variables
        self.p_otu=b2.shape[0]
        self.p_met=b3.shape[0]


    def simulate_mediation(self,x_func,m_func,x_args=None,m_args=None):
        """
        @param x_func: function to generate X (exposure)
        @param m_func: function to generate M (mediator)
        @param x_args: optional parameters to pass in x_func
        @param m_args: optional parameters to pass in m_func
        """
        x_shape=(self.p_otu,self.n)
        m_shape=(self.p_met,self.n)

        if x_args!=None:
            self.exposure=x_func(*x_shape,**x_args)
        else:
            self.exposure=x_func(*x_shape)
        
        if m_args!=None:
            self.mediator=m_func(*m_shape,**m_args)
        else:
            self.mediator=m_func(*m_shape)
        
        # simulate mediation effect
        
        self.mediator+=np.dot(self.b2.T,self.exposure)+np.random.normal(size=self.mediator.shape)
        self.outcome=np.dot(self.b3.T,self.mediator)+np.dot(self.b1.T,self.exposure)

        #add noise again
        self.outcome=self.outcome+np.random.normal(size=self.outcome.shape)


    def get_coeff(self):
        return {"b1":self.b1,"b2":self.b2,"b3":self.b3}


    def set_coeff(self,b1,b2,b3):
        self.b1=b1
        self.b2=b2
        self.b3=b3

    def estimate(self,solver):
        """
        estimate the parameters with B&K steps
        @param solver: solver to use for solving the linear systems
        """
        X=self.exposure.copy()
        Y=self.outcome.copy()
        M=self.mediator.copy()

        # step 1.
        X_1=np.vstack([X,np.ones(X.shape[1])])
        A1=solver(X_1.T,Y.T)
        a11=A1[:-1]

        # step 2.
        A2=solver(X_1.T,M.T)
        a21=A2[:-1]

        # step 3.
        X_M_1=np.vstack([X,M,np.ones(X.shape[1])])
        A3=solver(X_M_1.T,Y.T)
        a31=A3[:X.shape[0]]
        a32=A3[X.shape[0]:-1]

        self.A={'a11':a11,
                'a21':a21,
                'a31':a31,
                'a32':a32}

    def score(self):
        """
        calculates MSE of the prediction
        """
        self.A_true={'a11':self.b2.dot(self.b3)+self.b1,
                    'a21':self.b2,
                    'a31':self.b1,
                    'a32':self.b3}
        self.MSE=dict()

        for a in self.A_true.keys():
            mse=self.__mse(self.A[a],self.A_true[a])
            self.MSE[a]=mse

    def _truth(self):
        return {'a11':self.b2.dot(self.b3)+self.b1,
                    'a21':self.b2,
                    'a31':self.b1,
                    'a32':self.b3}


    def __mse(self,a,b):
        mse=np.power(b-a,2).mean()

        return mse

    def __joint_sig_test(self,a,b):
        #! what is the estimate of standard error for b?
        #! using rmse instead
        dist=norm(loc=0,scale=1)
        rmse=np.sqrt(np.power(b-a,2).mean())
        phi=np.abs(a)/rmse

        return 2*(1-dist.cdf(phi))

    def __pairwise_pr(self,X):
        corr=dict()
        ps=dict()
        shape=X.shape[0]

        for i in range(shape):
            for j in range(shape):
                if i!=j:
                    c,p=pearsonr(X[i],X[j])
                    corr[(i,j)]=c
                    ps[(i,j)]=p
        return corr,ps

    def pearsonCorr(self):
        """
        computes pairwise pearson correlation
        """
        
        return {
            "exposure":self.__pairwise_pr(self.exposure),
            "mediator":self.__pairwise_pr(self.mediator)
        }
