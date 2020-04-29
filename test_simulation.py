from numpy.linalg import inv,qr,solve,lstsq
from mediation_simulation import *

var={'n':2000,
     'p_otu':8, 
     'p_metabolite':4, 
     'mediations':1,
     'mediation_func':simple_mediation_direct
    }

truth_coef={
    'beta11':-2,
    'beta21':2,
    'beta32':-1.5
}

def generate(n=2000):
    var['n']=n
    output=simulate_simple_mediation(**var)
    producer,target=output["OTU"][:-1,:],output["OTU"][-1,:]
    metabolite=output["Metabolite"]
    
    # step 1
    X=np.vstack([np.ones((1,producer.shape[1])),producer])
    Y=target

    coef1,_,_,_=lstsq(X.T,Y.T,rcond=None)
    b10,b11=coef1[0],coef1[1:]
    
    
    # step 2
    Me=metabolite
    coef2,_,_,_=lstsq(X.T,Me.T,rcond=None)

    b20,b21=coef2[0],coef2[1:]
    
    
    # step 3
    X_Me=np.vstack([X,Me])
    coef3,_,_,_=lstsq(X_Me.T,Y.T,rcond=None)

    b30,b31,b32=coef3[0],coef3[1:X.shape[0]],coef3[X.shape[0]:]
    
    # results
    fitted = {
        'beta11':np.asarray(b11).reshape(-1),
        'beta21':np.asarray(b21),
        'beta31':np.asarray(b31).reshape(-1),
        'beta32':np.asarray(b32).reshape(-1)
    }
    
    temp=np.zeros((var["p_otu"]-1,var["p_metabolite"]))
    temp[-1,-1]=truth_coef["beta21"]
    
    true_vals={
        'beta11':np.hstack([np.zeros(X.shape[0]-2),truth_coef["beta32"]*truth_coef["beta21"]]),
        'beta21':temp,
        'beta31':np.hstack([np.zeros(X.shape[0]-1)]),
        'beta32':np.hstack([np.zeros(var["p_metabolite"]-var["mediations"]),truth_coef["beta32"]])
    }
    
    # sum of sqaured differences
    d={}
    for k in true_vals.keys():
        d[k]=np.power(fitted[k]-true_vals[k],2).sum()
    
    return fitted,true_vals,d

