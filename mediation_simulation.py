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

def multi_var_mediation(n, b11 = -2, b21 = 2, b32 = -1.5):
    """
    A is a symmetric positive semidefinite matrix,
    encoding relationships inbetween each pair of OTUs
    """
    A=_A(n)
    producer=(np.random.multivariate_normal(np.zeros(n),A)
             +np.random.rand(n))
    
    mm = b21*producer + np.random.normal(size=n)
    target = b11*producer + b32*mm + np.random.normal(size=n)

    otu = np.vstack([producer,target])
    temp = {"OTU": otu,"Metabolite": np.expand_dims(mm, axis=0)}

    return temp
    

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


