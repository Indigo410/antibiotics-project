import numpy as np
from numpy.linalg import qr,solve,inv,cholesky

def linmodEst(x,y):
    """
    reimplementation of linmodEst function from R
    """
    q,r=qr(x)
    xTx_1=inv(r.T.dot(r))
    coef=xTx_1.dot(x.T).T.dot(y.T).T

    # degrees of freedom & std of residuals
    df=x.shape[0]-x.shape[1]
    sigma2=np.power(y-coef.dot(x),2).sum()/df

    # sigma^2*(x'x)^-1
    vcov = sigma2*(xTx_1)

    return coef, vcov, np.sqrt(sigma2),df


