import numpy as np
from scipy.linalg import eigh


def _A(n,pairs=((1,0),(3,2))):
    """
    helper method for producing the cov matrix A
    for multivariate mediation simulation
    """
    A=np.zeros((n,n))
    np.fill_diagonal(A,1)
    for i,j in pairs:
        A[i,j]=1
        A[j,i]=1
    return A
    
def simple_mediation(n, beta11 = -2, beta21 = 2, b32 = -1.5):
    producer = np.random.rand(n)
    mm = beta21*producer + np.random.normal(size=n)
    target = beta11*producer + b32*mm + np.random.normal(size=n)

    otu = np.vstack([producer,target])
    temp = {"OTU": otu,"Metabolite": np.expand_dims(mm, axis=0)}

    return temp
    
def simple_mediation_direct(n, beta11 = -2, beta21 = 2, b32 = -1.5):
    producer = np.random.rand(n)
    mm = beta21*producer + np.random.normal(size=n)
    target = b32*mm + np.random.normal(size=n)

    otu = np.vstack([producer,target])
    temp = {"OTU":otu,"Metabolite":np.expand_dims(mm, axis=0)}

    return temp

def simulate_simple_mediation(n, p_otu, p_metabolite, mediations = 1,mediation_func=simple_mediation_direct):
    not_mediation_otu = p_otu-2*mediations
    not_mediation_metabolite = p_metabolite-mediations

    otutable = np.random.rand(not_mediation_otu,n)
    metabolitetable = np.random.rand(not_mediation_metabolite,n)

    for _ in range(mediations):
        temp = mediation_func(n)
        otutable = np.vstack([otutable,temp["OTU"]])
        metabolitetable = np.vstack([metabolitetable,temp["Metabolite"]])

    final = {
        "OTU": np.matrix(otutable),
        "Metabolite": np.matrix(metabolitetable) ,
        "OtuAnn": {"Species":["OTU{}".format(i) for i in range(1,p_otu+1)]},
        "MetAnn": {"Species":["M{}".format(j) if j<not_mediation_metabolite else "M{} Mediator".format(j) for j in range(1,p_metabolite+1)]}
        }

    return final


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
        self.n=n
        self.p_otu=p_otu,
        self.p_metabolite=p_metabolite
        self.mediations=mediations
        self.b=b
        self.mediation_func=mediation_func
        self.multivar=multivar
        if multivar:
            self.A=A
            self.B=B
            
        # ? unsure about the metabolite abundance
        if multivar:
            self.otutable = np.random.multivariate_normal(np.zeros(p_otu),A,size=(p_otu,n))
            self.metabolitetable = np.random.multivariate_normal(np.zeros(p_metabolite),B,size=(p_metabolite,n))
        else:
            self.otutable = np.random.rand(p_otu,n)
            self.metabolitetable = np.random.rand(p_metabolite,n)
    
    def simulate_abundance(self,
        mediations=1,
        b={'b11':-2,'b21':2,'b32':-1.5},
        mediation_func=self._simple)
        """
        mediations: int, number of mediations specified for the simulation
        b: dictionary, truth values for the relationship between producer and target
        mediation_func: callable, the method for generating related producer-target abundance
        """
        for i in range(mediations):
            otu,metabolite=mediation_func(self.p_otu-1-i,**self.b)
            self.otutable[ix,:]=otu
            self.metabolitetable[ix,:]=metabolite
        
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
