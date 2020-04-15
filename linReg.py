import numpy as np
from numpy.linalg import qr,solve,inv,cholesky

def linmodEst(x,y):
    """
    reimplementation of linmodEst function from R
    """
    q,r=qr(x)
    xTx_1=inv(r.T.dot(r))
    coef=xTx_1.dot(x.T).dot(y)

    # degrees of freedom & std of residuals
    df=x.shape[0]-x.shape[1]
    sigma2=sum((y-x*coef)**2)/df

    # sigma^2*(x'x)^-1
    vcov = sigma2*(xTx_1)

    return coef, vcov, np.sqrt(sigma2),df


