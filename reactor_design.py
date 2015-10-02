"""Contains reactor design functions"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def volume_CSTR(FA0, rate, Xmin, Xmax):
    volume = (FA0/-rate)*(Xmax-Xmin)
    return volume

def volume_PFR(FA0, rate, Xmin, Xmax):
    volume = compute_integral(np.linspace(Xmin,Xmax,np.size(rate)), FA0/-rate)
    return volume

def compute_integral(x_vals, y_vals):
    a = x_vals[:-1]
    b = x_vals[1:]    
    integral = 0.5*(b-a)*(y_vals[1:] + y_vals[:-1])
    integral = sum(integral)
    return integral

def plot_levenspiel(X, FA0_rA):
    
    plt.figure()
    plt.plot(X, FA0_rA)
    plt.xlabel("Conversion, X")
    plt.ylabel("FA0/-rA")
    plt.title("Levenspiel Plot")
    
    plt.savefig("levenspiel_plot.png")
    
    plt.close()