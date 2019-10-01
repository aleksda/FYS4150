import os	# Change dir to previous directory
os.chdir('/home/aleksandar/Desktop/FYS4150/Project2/')

import matplotlib.pyplot as plt
import pandas as pd, numpy as np
import seaborn as sns
sns.set()

#bc = pd.read_fwf('data/buckling.txt', header = None)
q1 = pd.read_fwf('data/quantumone.txt', header = None)
q2w01 = pd.read_fwf('data/quantumtwo_w01.txt', header = None)
q2w1 = pd.read_fwf('data/quantumtwo_w1.txt', header = None)
q2w5 = pd.read_fwf('data/quantumtwo_w5.txt', header = None)

#print(df)

#plt.plot(abs(bc.ix[:,0] ** 2))

plt.plot(abs(q1.loc[:,0] ** 2), label=r'$\psi_{v_{0}}$')
plt.plot(abs(q1.loc[:,1] ** 2), alpha=0.8, label=r'$\psi_{v_{1}}$')
plt.plot(abs(q1.loc[:,2] ** 2), alpha=0.8, label=r'$\psi_{v_{2}}$')

plt.xlabel(r'$\rho$')
plt.ylabel(r'radial solution $|u(\rho)|^{2}$')
plt.legend(loc='best')
plt.savefig('figures/quantumdot_one.jpg')
plt.show()

plt.plot(abs(q2w01.loc[:,0] ** 2), label=r'$\psi_{v_{0}}, \omega=0.1$')
plt.plot(abs(q2w1.loc[:,0] ** 2), label=r'$\psi_{v_{0}}, \omega=1.0$')
plt.plot(abs(q2w5.loc[:,0] ** 2), label=r'$\psi_{v_{0}}, \omega=5.0$')
plt.xlabel(r'$\rho$')
plt.ylabel(r'radial solution $|\psi(\rho)|^{2}$')
plt.legend(loc='best')
plt.savefig('figure/quantumdot_two.jpg')
plt.show()
