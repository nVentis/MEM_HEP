import numpy as np
from scipy.stats.qmc import Sobol
from scipy.integrate import quad
from matplotlib import pyplot as plt

n_iter = np.geomspace(10, 10000, dtype=int)

# Integrate f(x)=5x^2 between 3 and 5
# F(x) = 5/3x^3 + c
# I = 5/3*(125-27) = 5/3*(98) = 160

def f(x):
    #return 5 * x**2
    return 1/x

def mc_integrate(func, lower, upper, iterations):
    x = lower + np.random.rand(iterations) * (upper-lower)
    y = func(x)
    
    V = upper - lower
    
    return y.sum()*1/(iter) * V

def qmc_integrate(func, lower, upper, iterations):
    #print(iterations)
    rand = Sobol(1)
    iter_pow_of_2 = int(np.ceil(np.log(iterations)/np.log(2)))

    x = lower + rand.random_base2(iter_pow_of_2) * (upper-lower)
    # print("Pre:{}|Post:{}".format(len(x), iterations))
    x = x[0:iterations]
    
    y = func(x)
    
    V = upper - lower
    
    return y.sum()*1/(iter) * V

convs = []
qmcs = []

cum_convs = []
cum_qmcs = []

time_conv = []
time_qmc = []

bounds = (3,5)
corr, err = quad(f, bounds[0], bounds[1])

for i in range(len(n_iter)):
    iter = n_iter[i]
    
    result = mc_integrate(f, bounds[0], bounds[1], iter)
    convs.append(result)
    
    result = qmc_integrate(f, bounds[0], bounds[1], iter)
    qmcs.append(result)
    
    cum_conv = 0
    cum_qmc = 0
    n_cum = n_iter[:i].sum()
    
    # Weight by number of iterations
    for j in range(i):
        cum_conv += convs[j]*(n_iter[j]/n_cum)
        cum_qmc += qmcs[j]*(n_iter[j]/n_cum)
        
    cum_convs.append(cum_conv)
    cum_qmcs.append(cum_qmc)
        

print(corr)

show_cumulated = True

plt.plot(n_iter, cum_convs if show_cumulated else convs, label="MC")
plt.plot(n_iter, cum_qmcs if show_cumulated else qmcs, label="QMC")
plt.plot(n_iter, corr*np.ones(len(n_iter)), label="Correct")
plt.legend(loc="upper left")
plt.suptitle("Cumulated runs" if show_cumulated else "Individual runs")
plt.xscale("log")
plt.yscale("log")
plt.show()
