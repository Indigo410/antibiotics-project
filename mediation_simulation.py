import numpy as np

def simulate_simple_mediation(n, p_otu, p_metabolite, mediations = 1):
    not_mediation_otu = p_otu-2*mediations
    not_mediation_metabolite = p_metabolite-mediations

    otutable = np.random.rand(not_mediation_otu,n)
    metabolitetable = np.random.rand(not_mediation_metabolite,n)

    for _ in range(mediations):
        temp = simple_mediation_direct(n)
        otutable = np.vstack([otutable,temp["OTU"]])
        metabolitetable = np.vstack([metabolitetable,temp["Metabolite"]])

    final = {
        "OTU": np.matrix(otutable),
        "Metabolite": np.matrix(metabolitetable) ,
        "OtuAnn": {"Species":["OTU{}".format(i) for i in range(1,p_otu+1)]},
        "MetAnn": {"Species":["M{}".format(j) if j<not_mediation_metabolite else "M{} Mediator".format(j) for j in range(1,p_metabolite+1)]}
        }

    return final

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
