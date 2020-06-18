import numpy as np
from scipy.stats import norm,pearsonr
from Simulate import Simulate


def score(s,metric=_MSE,args={}) -> dict:
    """
    
    evaluates the simulator with given metric,
    metric is default to MSE
    """
    if args:
        result=metric(s,**args)
    else:
        result=metric(s)
    
    return result
    

def _MSE(s) -> float:
    """
    computes the MSE value between truth vector and estimated vector
    """
    truth=s._truth()
    mse=dict()

    for k in truth.keys():
        mse[k]=np.power(truth[k],s.A[k],2).mean()

    return mse



################################################
#               metric fucntions               #
# input : simulator, (additonal args,)         #
# output: p values for all estimated vectors   #
################################################


def non_param_bootstrap(s,solver,n=2000,level=0.95) -> dict:
    """
    computes p-value of the estimation by non-parametric bootstrap
    @param n: number of trials for bootstrapping,default to 2000
    @param level: confidence level, default to 0.95
    """
    null_dist=dict=(a11=np.zeros((s.p_otu,n)),
                    a21=np.zeros((s.p_met,n)),
                    a31=np.zeros((s.p_otu,n)),
                    a32=np.zeros((s.p_met,n)))
    std_dist=null_dist.copy()

    # begin bootstrapping
    for i in range(n):
        # null 
        # renomalization with randomized outcome abundance
        outcome=np.random.rand(s.n)
        normalized=np.vstack([s.exposure,outcome])
        normalized/=np.sum(normalized,axis=0)
        normalized*=np.vstack([s.exposure,s.outcome]).sum()

        exposure_n=normalized[:s.exposure.shape[0]]
        outcome_n=normalized[s.exposure.shape[0]:]

        # solve for coefficients
        A_null=s.B_K_steps(solver,X=exposure_n,Y=outcome_n,M=s.mediator)


        # standard bootstrap
        # shuffle horizontally accross each trial
        rng=np.random.default_rng()
        exposure=s.exposure.copy()
        rng.shuffle(exposure,axis=1)

        A_std=s.B_K_steps(solver,X=exposure,Y=s.outcome,M=s.mediator)

        # store current iteration result
        for k in null_dist.keys():
            null_dist[k][i]=A_null[k]
            std_dist[k][i]=A_std[k]

    # TODO: 2 sample t test



    raise NotImplementedError


################################################
#          other correlation metrics           #
# input : simulator, (additonal args,)         #
# output: correlation coefficients and         #
#         p values for OTU and mediators       #
################################################

def pearson_corr(simulator,threshold=0.01):
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

        
def joint_significance(simulator):
    raise NotImplementedError
