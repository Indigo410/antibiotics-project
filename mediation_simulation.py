import numpy as np
import pandas as pd

def simulate_simple_mediation(n, p_otu, p_metabolite, mediations = 1):
    not_mediation_otu = p_otu-2*mediations
    not_mediation_metabolite = p_metabolite-mediations

    otutable = np.random.rand(not_mediation_otu,n)
    metabplitetable = np.random.rand(not_mediation_metabolite,n)

    # format the otutable for concatenation
    otutable = pd.DataFrame(otutable,columns=("producer","target"))

    for m in range(1, mediations):
        temp = simple_mediation_direct(n)
        otutable = pd.concat([otubale,temp["OTU"]], "inner")  # inner join for mimicing the behaviour of rbind
        p_metabolite = np.concatenate([metabplitetable,temp["Metabolite"]])

    final = {
        "OTU": otutable,
        "Metabolite": metabplitetable,
        "OtuAnn": pd.DataFrame({"Species":["OTU{}".format(i) for i in range(1,p_otu+1)]}),
        "MetAnn": pd.DataFrame({"Species":["M{}".format(j) if j<not_mediation_metabolite else "M{} Mediator".format(j) for j in range(1,p_metabolite+1)]})
        }

    return final

def simple_mediation(n, beta11 = -2, beta21 = 2, b32 = -1.5):
    producer = np.random.rand(n)
    mm = beta21*producer + np.random.normal(size=n)
    target = beta11*producer + b32*mm + np.random.normal(size=n)

    # format the results for output
    otu = pd.DataFrame(data={"producer":producer,"target":target})
    temp = {"OTU": otu,"Metabolite": mm}

    return temp
    
def simple_mediation_direct(n, beta11 = -2, beta21 = 2, b32 = -1.5):
    producer = np.random.rand(n)
    mm = beta21*producer + np.random.normal(size=n)
    target = b32*mm + np.random.normal(size=n)

    # format the results for output
    otu = pd.DataFrame(data={"producer":producer,"target":target})
    temp = {"OTU":otu,"Metabolite":mm}

    return temp
