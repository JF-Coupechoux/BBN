import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(18,13))

plt.subplot(121)
fichier= np.loadtxt('ML.out')
fichier4= np.loadtxt('MLL.out')
axes = plt.gca()
plt.grid(True)
plt.xlabel('Scale factor')
plt.ylabel('Evolution of energy density')
plt.plot(fichier[:,1],fichier[:,5], label=r'$\Omega_{b}$ baryon density')
plt.plot(fichier[:,1],fichier[:,4], label=r'$\Omega_{r}$ radiation density')
plt.plot(fichier[:,1],fichier[:,7], label=r'$\Omega_{\Lambda}$ dark energy density')
plt.plot(fichier4[:,1],fichier4[:,6], label=r'$\Omega_{\phi}$ field density $(m=1.5\times 10^{-20}eV$, $\lambda=10^{-81})$')
plt.yscale('log')
plt.xscale('log')
axes.set_xlim(1.e-12,1)
axes.set_ylim(1.e-4,1.e44)
plt.legend(loc='upper right')

plt.subplot(122)
fichier= np.loadtxt('ML.out')
fichier1= np.loadtxt('ML1.out')
fichier2= np.loadtxt('ML2.out')
axes = plt.gca()
plt.grid(True)
plt.xlabel('Scale factor')
plt.ylabel('Evolution of energy density')
plt.plot(fichier[:,1],fichier[:,5], label=r'$\Omega_{b}$ baryon density')
plt.plot(fichier[:,1],fichier[:,4], label=r'$\Omega_{r}$ radiation density')
plt.plot(fichier[:,1],fichier[:,7], label=r'$\Omega_{\Lambda}$ dark energy density')
plt.plot(fichier[:,1],fichier[:,6], label=r'$\Omega_{\phi}$ field density $(m=10^{-24}eV$, $\lambda=10^{-98})$')
plt.plot(fichier1[:,1],fichier1[:,6], label=r'$\Omega_{\phi}$ field density $(m=10^{-24}eV$, $\lambda=10^{-100})$')
plt.plot(fichier2[:,1],fichier2[:,6], label=r'$\Omega_{\phi}$ field density $(m=10^{-24}eV$, $\lambda=10^{-111})$')
plt.yscale('log')
plt.xscale('log')
axes.set_xlim(1.e-12,1)
axes.set_ylim(1.e-4,1.e44)
plt.legend(loc='upper right')

plt.show()
