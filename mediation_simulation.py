import numpy as np
from scipy.linalg import eigh


def multi_var_mediation(n, b11 = -2, b21 = 2, b32 = -1.5, A=None):
    # ? not sure how this psd covariance matrix is supposed to be
    if A==None:
        psd=False
        while not psd:
            # randomly generate a psd matrix
            mat=np.random.rand(n,n)
            mat=mat.T.dot(mat)

            # check whether mat is indeed psd
            E,_=eigh(mat)
            psd=all(E>=0)
    else:
        mat=A.T.dot(A)

    producer=(np.random.multivariate_normal(np.zeros(n),mat,size=n)
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


