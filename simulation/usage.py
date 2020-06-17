from Simulate import Simulate,coefficient
from Random import TreeSimilarity
from Test import *
import pycasso
import numpy as np
import pandas as pd    # for visualizing result

# initialize variables
n=500
p_otu=3
p_met=10
b1=np.zeros((p_otu,1))
b2=np.zeros((p_otu,p_met))
b3=np.zeros((p_met,1))

b1[0,-1]=-2
b2[0,-1]=2
b3[0,-1]=-1.5

b=dict(b1=b1,b2=b2,b3=b3)
var=[n,b]


# function declaration
def lasso(x,y):
    sol = pycasso.Solver(x,y)
    sol.train()
    return sol.coef()["beta"][-1]

def pipeline(random_state=0):
    if random_state is None:
        pass
    else:
        np.random.seed(random_state)
    
    s=Simulate(*var)
    t_exposure=TreeSimilarity(p_otu)
    t_mediator=TreeSimilarity(p_met)
    s.simulate_mediation(t_exposure.rand,t_mediator.rand)
    s.estimate(lasso)
    
    evaluation=score(s,non_param_bootstrap)

    return evaluation


if __name__=="__main__":
    evaluation=pipeline()
    
    for k in evaluation.keys():
        print(k,evaluation[k])