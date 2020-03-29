import numpy as np

def simulate_simple_mediation(n, p_otu, p_metabolite, mediations = 1):
    not_mediation_otu = p_otu-2*mediations
    not_mediation_metabolite = p_metabolite-mediations
    otutable = np.random.rand(not_mediation_metabolite,n)
    metabplitetable = np.random.rand(not_mediation_metabolite,n)

    for m in range(1, mediations):
        # TODO: fix this after finishing simple_mediation_direct
        # NOTE: temp is not a data frame
        temp = simple_mediation_direct(n)  # assuming temp is a data frame
        otutable = pd.concat([otubale,temp['OTU']], "inner")  # inner join for mimicing the behaviour of rbind
        p_metabolite = pd.concat([metabolitetable, temp['Metabolite']], "inner")
    
    # TODO

def simple_mediation(n, beta11 = -2, beta21 = 2, b32 = -1.5):
    producer = np.random.rand(n)
    mm = beta21*producer + np.random.normal(size=n)
    target = beta11*producer + b32*mm + np.random.normal(size=n)

    # format the results for output
    otu = np.concatenate((producer, target))
    temp = {"OTU": otu,"Metabolite": mm}

    return temp
    
def simple_mediation_direct(n, beta11 = -2, beta21 = 2, b32 = -1.5):
    producer = np.random.rand(n)
    mm = beta21*producer + np.random.normal(size=n)
    target = b32*mm + np.random.normal(size=n)
