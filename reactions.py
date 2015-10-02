"""General-ish stuff for modeling reaction kinetics. Currently only configured for liquid-phase power models."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

###REACTION CLASS HOPEFULLY TO BE USED AND EXPANDED AT SOME POINT###
class reaction_power_model(object):
    def __init__(self, rate_constant, powers):
        self.rate_constant = rate_constant
        self.reactant_orders = powers
        self.overall_order = sum(powers)

    def power_rate_law(self, conc):
        k = self.rate_constant
        powers = self.reactant_orders
        rate = k*np.prod(np.power(conc, powers))
        return rate

    def display_info(self):
        print "Rate constant is: ", self.rate_constant
        print "Reactant orders are: ", self.reactant_orders
        print "Overall reaction order is: ", self.overall_order

###Functions that should (hopefully) function for now###
def fit_power_model(conc, rates):
    # Use linear least squares to fit rate vs concentration data
    # Column vectors of conc should correspond to concentrations of each reactant
    
    npoints = np.size(rates)
    
    rates = np.reshape(rates, (npoints,1))
    
    if np.size(np.shape(conc)) == 1:
        conc = np.reshape(conc, (npoints,1))
    elif not np.shape(conc)[0] == npoints:
        conc = np.reshape(conc, (npoints, np.shape(conc)[0]))
        
    X = np.ones((np.shape(conc)[0], np.shape(conc)[1]+1))
    X[:,1:] = conc
    print np.shape(X)
    
    M = np.dot(np.transpose(X), X)
    b = np.dot(np.transpose(X), rates)
    
    params = np.linalg.solve(M, b)
    
    k = np.exp(params[0])
    powers = params[1:]
    
    return params, k, powers

def get_rates_fd(conc, step):
    # Use finite difference to approximate rates
    
    nstep = np.size(conc)
    conc = np.reshape(conc, (nstep,))
    weights = np.array((-3,4,-1))
    rates = np.zeros(np.shape(conc))
    rates[0] = np.dot(weights, conc[:3])/(2*step)
    rates[-1] = np.dot(np.flipud(-weights), conc[-3:])
    
    rates[1:-1] = (conc[2:] - conc[0:-2])/(2*step)
    
    return rates

def power_model(t, conc, k, powers):
    #Same as power_rate_law from reaction_power_model class
    rate = -k*np.prod(np.power(conc, powers))
    return rate

def integrate_rate_law(t_range, step, initial_cond, k, powers):
    # Numerically integrate rate law using 4th order Runge-Kutta
    n_reactants = np.size(initial_cond)
    x = np.tile(t_range, (n_reactants, 1))
    y = np.zeros((n_reactants, np.shape(t_range)[0]))
    y[:,0] = initial_cond
    
    weights = np.transpose(np.array([1,1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0], dtype=float))
    
    for n in range(np.size(t_range)-1):
        k1 = step*power_model(x[:,n], y[:,n], k, powers)
        k2 = step*power_model(x[:,n]+0.5*step, y[:,n]+0.5*k1, k, powers)
        k3 = step*power_model(x[:,n]+0.5*step, y[:,n]+0.5*k2, k, powers)
        k4 = step*power_model(x[:,n]+step, y[:,n]+k3, k, powers)
        
        y[:,n+1] = np.dot(np.array((y[:,n], k1, k2, k3, k4)), weights)

    return y

def plot_time(t_range, concentration, xlabel, ylabel, title, filename):
    # Plots concentration vs time
    
    plt.figure()
    plt.plot(t_range, concentration)

    plt.xlabel("%s"%xlabel)
    plt.ylabel("%s"%ylabel)
    plt.title("%s"%title)

    plt.savefig("%s.png"%filename)
    plt.close()