import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor,DistanceMatrix

class TreeSimilarity:
    """
    generate covariance matrix as a parameter of multinomial distribution
    from a tree of similarities between OTUs
    """
    
    def __init__(self,n,names=None):
        if names is None:
            names=["otu{}".format(i) for i in range(n)]
        else:
            names=[names+str(i) for i in range(n)]
        self.names=names
        
        self.distMat=np.tril(np.random.rand(n,n))

        temp=np.tril(self.distMat,-1)
        # consider covariance as the inverse of distance
        self.cov=1-(temp+temp.T)
        self.cov=np.cov(self.cov)

    def draw(self):
        """
        visualize the phylo tree
        """
        mat=list(map(lambda x: list(filter(lambda x:x>0,x)),self.distMat.tolist()))
        constructor = DistanceTreeConstructor()
        upgmatree = constructor.upgma(DistanceMatrix(self.names,mat))
        
        Phylo.draw_ascii(upgmatree)

    def rand(self,sizeX,sizeY):
        return np.random.multivariate_normal(np.zeros(sizeX),self.cov,size=sizeY).T