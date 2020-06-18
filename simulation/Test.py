import numpy as np
from scipy.stats import norm,pearsonr
from Simulate import Simulate

def _MSE(s) -> dict:
    """
    computes the MSE value between truth vector and estimated vector
    """
    truth=s._truth()
    mse=dict()

    for k in truth.keys():
        mse[k]=np.power(truth[k],s.A[k],2).mean()

    return mse

def score(simulator,metric=_MSE,args={}) -> dict:
    """
    general function for evaluating the estimation.
    evaluates the simulator with given metric,
    @param simulator: the Simulate object to be evaluated
    @param metric: metric to evaluate, default to MSE
    @param args: additional arguments to pass into metric function
    """
    if args:
        result=metric(simulator,**args)
    else:
        result=metric(simulator)
    
    return result
    

################################################
#               metric fucntions               #
# input : simulator, (additonal args,)         #
# output: p values for all estimated vectors   #
################################################


def non_param_bootstrap(s,solver,n=1000) -> dict:
    """
    computes p-value of the estimation by non-parametric bootstrap
    @param n: number of trials for bootstrapping,default to 2000
    """
    null_dist=dict(a11=np.zeros((s.p_otu,n)),
                    a21=np.zeros((s.p_otu,n)),
                    a31=np.zeros((s.p_otu,n)),
                    a32=np.zeros((s.p_met,n)))
    std_dist=dict(a11=np.zeros((s.p_otu,n)),
                    a21=np.zeros((s.p_otu,n)),
                    a31=np.zeros((s.p_otu,n)),
                    a32=np.zeros((s.p_met,n)))

    # begin bootstrapping
    for i in range(n):
        #----null----#
        # renomalization with randomized outcome abundance
        outcome=np.random.rand(s.n)
        normalized=np.vstack([s.exposure.copy(),outcome])
        normalized/=np.sum(normalized,axis=0)
        normalized*=np.vstack([s.exposure.copy(),s.outcome.copy()]).sum()

        exposure_n=normalized[:s.p_otu]
        outcome_n=normalized[s.p_otu:]

        # solve for coefficients
        A_null=s.B_K_steps(solver,X=exposure_n,Y=outcome_n,M=s.mediator)

        #----standard bootstrap----#
        rng = np.random.default_rng()
        ix=rng.choice(s.n,s.n)

        A_std=s.B_K_steps(solver,X=s.exposure[:,ix],Y=s.outcome[:,ix],M=s.mediator)

        # store current iteration result
        for k in A_null.keys():
            null_dist[k][:,i]=A_null[k]
            std_dist[k][:,i]=A_std[k]


    result=dict()
    for k in null_dist.keys():
        cnt=((null_dist[k]-std_dist[k])!=0).sum(axis=1)
        result[k]=cnt/n       
        
    return result



################################################
#          other correlation metrics           #
# input : simulator, (additonal args,)         #
# output: correlation coefficients and         #
#         p values for OTU and mediators       #
################################################

def pearson_corr(simulator,threshold=0.01)->dict:
    """
    pairwise pearson correlation between all OTU and mediators
    """
    result=dict(exposure=dict(),mediator=dict())

    for k in ("exposure","mediator"):
        corr=dict()
        ps=dict()
        X=getattr(simulator,k)
        shape=X.shape[0]

        for i in range(shape):
            for j in range(shape):
                if i!=j:
                    c,p=pearsonr(X[i],X[j])
                    
                    # keep only the ones with statistic significance
                    if p>=threshold:
                        corr[(i,j)]=c
                        ps[(i,j)]=p

        result[k]=dict(corr=corr,p=ps)
    
    return result

        
def joint_significance(simulator):
    raise NotImplementedError
