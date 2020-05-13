import numpy as np


class Simulation:
    def __init__(self,
        n,
        p_otu,
        p_metabolite,
        multivar=False,
        A=None,
        B=None):
        """
        n: int, number of simulation samples
        p_otu: int, number of otu populations
        p_metabolite: int, number of metabolite populations
        multivar: bool, flag for whether the otu populations are correlated
        A: a (p_otu*p_otu) covariance matrix encoding the relationship between otu species; ignored
        if mutlivar is False; required if multivar is True
        B: a (p_metabolite*p_metabolite) covariance matrix
        """
        self.n=n
        self.p_otu=p_otu
        self.p_metabolite=p_metabolite
        self.multivar=multivar
        if multivar:
            self.A=A
            self.B=B
            
        if multivar:
            self.otutable = np.random.multivariate_normal(np.zeros(p_otu),A,size=(p_otu,n))
            self.metabolitetable = np.random.multivariate_normal(np.zeros(p_metabolite),B,size=(p_metabolite,n))
        else:
            self.otutable = np.random.rand(p_otu,n)
            self.metabolitetable = np.random.rand(p_metabolite,n)
    
    def run_simulation(self,
        mediations=1,
        b={'b11':-2,'b21':2,'b32':-1.5},
        mediation_func="simple"):
        """
        mediations: int, number of mediations specified for the simulation
        b: dictionary, truth values for the relationship between producer and target
        mediation_func: str or callable, the method for generating related producer-target abundance
        """
        self.b=b

        if isinstance(mediation_func,str):
            if mediation_func=="simple_direct": mediation_func=self.__simple_direct
            elif mediation_func=="simple": mediation_func=self.__simple
        elif callable(mediation_func):
            pass
        else:
            mediation_func=self.__simple

        self.is_dependent=np.zeros(self.otutable.shape[0])
        for ix in range(mediations):
            otu,metabolite=mediation_func(self.p_otu-1-ix,**self.b)
            self.otutable[ix,:]=otu
            self.metabolitetable[ix,:]=metabolite

            self.is_dependent[ix]=1
        
        return self.otutable,self.metabolitetable
        
    def __simple(self,j,b11,b21,b32):
        """
        private method for calculating related target & metabolite abundance
        """
        producer=self.otutable[j,:]
        mm=b21*producer+np.random.normal(size=self.n)
        target=b11*producer+b32*mm+np.random.normal(size=self.n)

        return target,mm

    def __simple_direct(self,j,b11,b21,b32):
        """
        private method for calculating directly related target & metabolite abundance
        """
        producer=self.otutable[j,:]
        mm=b21*producer+np.random.normal(size=self.n)
        target=b32*mm+np.random.normal(size=self.n)

        return target,mm


class TreeSimilarity:
    """
    generate covariance matrix as a parameter of multinomial distribution
    from a tree of similarities between OTUs
    """
    # TODO
    
    def __init__(self):
        pass

    def generate(self):
        pass