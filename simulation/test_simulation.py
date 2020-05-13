from numpy.linalg import inv,qr,solve,lstsq
from simulation import simulation

def generate(vars,n=2000):
    """
    generate n independent simulation trials, n default to 2000
    @param n: number of simulations
    @param vars: variables to pass into Simulation
    """
    pass
    # TODO

def truth_value(sim):
    """
    calculates the truth value for the given Simulation object
    @param sim: a Simulation object
    """
    pass
    # TODO

def b_k_steps(sim,truth):
    """
    calculate estimators of significant mediation effects
    and corresponding p-values using b&k steps
    @param sim: a Simulation object
    @param truth: corresponding truth values
    """
    pass
    # TODO
