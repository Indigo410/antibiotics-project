import numpy as np
from scipy.stats import norm
from numpy.linalg import inv,qr,solve,lstsq
from simulation import Simulation

def generate(vars0,vars1,t=2000):
    """
    generate t independent simulation trials
    @param t: number of simulations
    @param vars0: initializing variables to pass into Simulation
    @param vars1: varaibles to run the simulation
    """
    f={c:[] for c in ["beta11","beta21","beta31","beta32"]}

    for i in range(t):
        simulator=Simulation(**vars0)
        _=simulator.run_simulation(**vars1)
        truth=truth_value(simulator)
        
        fitted=b_k_steps(simulator,truth)
        for c in fitted.keys():
            f[c].append(fitted[c])
    return f

def run_test(differences):
    """
    run statistical test on the differences found
    for the coefficients predicted and truth values
    """
    pass
    #todo

def joint_significance(estimator,std):
    """
    joint significance test as defined in the paper btw351
    """
    dist=norm(loc=0,scale=1)
    phi=np.linalg.norm(estimator)/std

    return 2*(1-norm.cdf(phi))

def truth_value(sim):
    """
    calculates the truth value for the given Simulation object
    @param sim: a Simulation object
    """
    beta=sim.b
    temp=np.zeros((sim.p_otu-sim.mediations,sim.p_metabolite))
    temp[-1,-1]=beta["b21"]

    true_coeff={
        'beta11':np.hstack([np.zeros(sim.p_otu-2),beta["b32"]*beta["b21"]]),
        'beta21':temp,
        'beta32':np.hstack([np.zeros(sim.p_metabolite-sim.mediations),
                            beta["b32"]])
    }

    return true_coeff

def b_k_steps(sim,truth):
    """
    calculate estimators of significant mediation effects
    and corresponding p-values using b&k steps
    @param sim: pre-run simulation object
    @param truth: truth vector
    """
    otutable=sim.otutable
    metabolitetable=sim.metabolitetable

    producer,target=otutable[~sim.is_dependent],otutable[sim.is_dependent]

    # step 1
    X=np.vstack([np.ones((1,producer.shape[1])),producer])
    Y=target

    coef1,_,_,_=lstsq(X.T,Y.T,rcond=None)
    b10,b11=coef1[0],coef1[1:]

    # step 2
    coef2,_,_,_=lstsq(X.T,metabolitetable.T,rcond=None)

    b20,b21=coef2[0],coef2[1:]

     # step 3
    X_Me=np.vstack([X,metabolitetable])
    coef3,_,_,_=lstsq(X_Me.T,Y.T,rcond=None)

    b30,b31,b32=coef3[0],coef3[1:X.shape[0]],coef3[X.shape[0]:]

    # results
    fitted = {
        'beta11':np.asarray(b11).reshape(-1),
        'beta21':np.asarray(b21),
        'beta32':np.asarray(b32).reshape(-1)
    }
    
    # ! change this to adapt to the statistical tests decided later
    # # absolute difference
    # d={}
    # for k in fitted.keys():
    #     d[k]=np.abs(truth[k]-fitted[k])

    return fitted
