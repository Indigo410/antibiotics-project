import numpy as np
from scipy.stats import norm,pearsonr


class Simulate:
    def __init__(self,n:int,b:dict,vars=None) -> None:
        """
        @param n: number of samples
        @param b1: all 0 for direct mediation, some otu directly influenced target otu otherwise
        """
        # input variables
        self.n=n
        self.b=b

        # inferred variables
        self.p_otu=b["b2"].shape[0]
        self.p_met=b["b3"].shape[0]

    def simulate_mediation(self,x_func,m_func,x_args=None,m_args=None) -> None:
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
        self.mediator+=np.dot(self.b["b2"].T,self.exposure)+np.random.normal(size=self.mediator.shape)
        self.outcome=np.dot(self.b["b3"].T,self.mediator)+np.dot(self.b["b1"].T,self.exposure)

        #add noise again
        self.outcome=self.outcome+np.random.normal(size=self.outcome.shape)

        # ? normalize all abundance
        # normalized=np.vstack([self.exposure,self.outcome])
        # normalized/=np.sum(normalized,axis=0)
        
        # self.exposure_n=normalized[:self.exposure.shape[0]]
        # self.outcome_n=normalized[self.exposure.shape[0]:]

    def get_coeff(self) -> dict:
        return self.b

    def set_coeff(self,b:dict) -> None:
        self.b=b

    def get_truth(self) -> dict:
        return dict(a11=self.b["b2"].dot(self.b["b3"])+self.b["b1"].reshape(-1),
                    a31=self.b["b1"].reshape(-1),
                    a21=self.b["b3"].reshape(-1),
                    a32=self.b["b3"].reshape(-1))


    def estimate(self,solver) -> None:
        """
        estimate the parameters with B&K steps
        @param solver: solver to use for solving the linear systems
        """
        self.A=Simulate.B_K_steps(solver,self.exposure,self.mediator,self.outcome)

    @staticmethod
    def __step1(solver,X,Y) ->np.ndarray:
        X_1=np.vstack([X,np.ones(X.shape[1])])
        A1=solver(X_1.T,Y.T)
        
        return A1[:-1]

    @staticmethod
    def __step2(solver,X,M) ->np.ndarray:
        X_1=np.vstack([X,np.ones(X.shape[1])])
        A2=solver(X_1.T,M.T)

        return A2[:-1]

    @staticmethod
    def __step3(X,M,Y) ->tuple:
        X_M_1=np.vstack([X,M,np.ones(X.shape[1])])
        A3=solver(X_M_1.T,Y.T)

        return A3[:X.shape[0]],A3[X.shape[0]:-1]
    
    @staticmethod
    def B_K_steps(solver,X:np.ndarray,M:np.ndarray,Y:np.ndarray) ->dict:
        """
        estimate the coefficients by B&K steps
        """
        a11=Simulate.__step1(solver,X,Y)

        a21=Simulate.__step2(solver,X,Y)

        a31,a32=Simulate.__step3(solver,X,M,Y)

        # format results
        A=dict(a11=a11,a21=a21,a31=a31,a32=a32)

        return A