""""
Script for plotting the power spectrum, variance and halo mass function in PBHs scenarios
Author: Pablo Villanueva Domingo
Started: October 2018
Last update: March 2021
"""
import matplotlib.pyplot as plt
from Source.functions import *

#--- POWER SPECTRUM PLOTS ---#

plot_pow = 0
plot_sigma = 0
plot_dndlnM = 1

cols = ["c","m","b"]
#cols = ["#B388EB","#8093F1","#9C0D38"]
MassPBH = 1.
#fracs = [1e2,1e1,1e0]
fracs = [1e-1,1e-2,1e-3]

z = 10

k = np.logspace(-3,4)
M = np.logspace(2,8,num=100)

fig_pow, (ax_pow) = plt.subplots(1,1)
fig_sig, (ax_sig) = plt.subplots(1,1)
fig_halo, (ax_halo) = plt.subplots(1,1)

if plot_pow==1:

    ax_pow.set_xlabel('$k \; [h/Mpc]$')
    #plt.set_ylabel('$P(k) \; [(Mpc/h)^3]$')
    ax_pow.set_ylabel('$\Delta^2(k)$')
    #plt.loglog(k, cosmo.matterPowerSpectrum(k), '-',label=r"$\Lambda CDM$")
    ax_pow.loglog(k, (k**3./(2.*np.pi**2.))*cosmo.matterPowerSpectrum(k)*cosmo.growthFactor(z)**2., 'k-',label=r"$\Lambda CDM$");
    for fpbh in fracs:
        #plt.loglog(k, PkPBH(k,fpbh,1.),color = cols[fracs.index(fpbh)], linestyle="-",label=r"$f_{PBH}M_{PBH}/M_{\odot}=$"+str(fpbh))
        ax_pow.loglog(k, (k**3./(2.*np.pi**2.))*PkPBH(k,fpbh,MassPBH)*cosmo.growthFactor(z)**2.,color = cols[fracs.index(fpbh)], linestyle="-",label=r"$f_{PBH}M_{PBH}/M_{\odot}=$"+scinot(fpbh))
    ax_pow.legend()
    ax_pow.set_xlim(1e-2,1e4)
    ax_pow.set_ylim(1e-4,1e5)
    ax_pow.set_title(r"$z=$"+str(z))
    fig_pow.savefig("Plots/powerspectrum.pdf", bbox_inches='tight')

if plot_sigma==1:

    ax_sig.set_xlabel('$M \; [M_{\odot}]$')
    ax_sig.set_ylabel('$\sigma(M)$')
    ax_sig.loglog(M,np.sqrt(sigvec(M,1e-10,MassPBH))*cosmo.growthFactor(z), 'k-',label=r"$\Lambda CDM$");
    for f, fpbh in enumerate(fracs):
        ax_sig.loglog(M,np.sqrt(sigvec(M,fpbh,MassPBH))*cosmo.growthFactor(z),color = cols[f], linestyle="-",label=r"$f_{PBH}M_{PBH}/M_{\odot}=$"+scinot(fpbh));
    ax_sig.legend(loc="upper right")
    ax_sig.set_xlim(M[0],M[-1])
    #plt.set_ylim(1e-10,1e6)
    ax_sig.set_title(r"$z=$"+str(z))
    fig_sig.savefig("Plots/sigma_pbh.pdf", bbox_inches='tight')

if plot_dndlnM==1:

    ax_halo.set_xlabel('$M \; [M_{sun}/h]$')
    ax_halo.set_ylabel('$dn/dlnM \; [(h/Mpc)^3]$')
    ax_halo.loglog(M, dndlnM_CDM(M,z), 'k-',label=r"$\Lambda CDM$")
    #np.savetxt("dndM_z_{:.1e}_LCDM_MJeans_{:.2e}_Mvir_{:.2e}.dat".format(z,MJeans(z,Tk_ad(z)),Mmin(1.e4,z)), np.transpose([M,dndlnM_CDM(M,z)]))
    for fpbh in fracs:
        ax_halo.loglog(M, dndlnM_PBHvec(M,z,fpbh,MassPBH), color = cols[fracs.index(fpbh)], linestyle="-",label=r"$f_{PBH}M_{PBH}/M_{\odot}=$"+scinot(fpbh))
        #np.savetxt("dndM_z_{:.1e}_MfPBH_{:.2e}_MJeans_{:.2e}_Mvir_{:.2e}.dat".format(z,fpbh,MJeans(z,Tk_ad(z)),Mmin(1.e4,z)), np.transpose([M,dndlnM_PBHvec(M,z,fpbh,MassPBH)]))
    ax_halo.axvline(x=MJeans(z,Tk_ad(z)),color="grey",linestyle=":")
    ax_halo.axvline(x=Mmin(1.e4,z),color="grey",linestyle="--")
    ax_halo.text(0.5*MJeans(z,Tk_ad(z)), 1e3, r"$M_{\rm J}(T_k = T_{ad})$", fontsize=10,rotation=90)
    ax_halo.text(0.5*Mmin(1.e4,z), 1e3, r"$M_{\rm vir}^{min}(T_{\rm vir} = 10^4 K)$", fontsize=10,rotation=90)
    ax_halo.legend(loc="lower left")
    ax_halo.set_xlim(M[0],M[-1])
    #plt.set_ylim(1e-10,1e6)
    ax_halo.set_title(r"$z=$"+str(z))
    fig_halo.savefig("Plots/halospbh.pdf", bbox_inches='tight')
