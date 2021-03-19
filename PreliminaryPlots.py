""""
Plot the gas density, spin temperature and optical depth as a function of radius
Pablo Villanueva Domingo
Started: October 2018
Last update: March 2021
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Source.functions import *


zvec = [10, 15]

xvec = np.logspace(-3,0,200)    # r/Rvir
#Ms = [1.e4,1.e5,1.e6,1.e7,1.e8] # Msun/h
Ms = [1.e4,1.e6,1.e8] # Msun/h
colors = ["c","m","b","g","r"]
lines = ["-",":"]


fig_rho, (ax_rho) = plt.subplots(1,1)
fig_Ts, (ax_Ts) = plt.subplots(1,1)
fig_tau, (ax_tau) = plt.subplots(1,1)

for iz, z in enumerate(zvec):

    for mm, M in enumerate(Ms):
        y = con(M,z)
        r = xvec*Rvir(M,z)
        Tk = Tvir(M,z)
        nH = rho_gas(M,z,y,r)*(1.-Y_He)/m_p
        Tspin = Ts(z,nH,Tk)

        logaloverRv = np.linspace(np.log(1.e-2),0.,num=200)    # al/Rvir(M,z)
        aloverRv = np.exp(logaloverRv)
        OpDepArr = []
        for a in aloverRv:
            OpDepArr.append(OptDepth21(M,z,y,a*Rvir(M,z)))

        #ax_rho.loglog(xvec, rho_nfw(M,z,y,r)/rho_c_z(z), linestyle=":", color=colors[mm])
        ax_rho.loglog(xvec, rho_gas(M,z,y,r)/rho_c_z(z), linestyle=lines[iz], color=colors[mm])
        ax_Ts.loglog(xvec, Tspin, linestyle=lines[iz], color=colors[mm])
        ax_Ts.loglog([xvec[0],xvec[-1]],[Tk,Tk],color=colors[mm],linestyle="--")
        ax_tau.loglog(aloverRv, OpDepArr, linestyle=lines[iz], color=colors[mm])


leg = []
for mm, M in enumerate(Ms):
    leg.append( Line2D([0], [0], color=colors[mm], linestyle="-", lw=4, label=r"$M_h=$"+scinot(M)+"$ M_{\odot}/h$") )
for iz, z in enumerate(zvec):
    leg.append( Line2D([0], [0], color="k", linestyle=lines[iz], lw=2, label="$z=$"+str(z) ) )

ax_rho.set_xlim(xvec[0],xvec[-1])
#ax_rho.set_ylim(1e-2,1e8)
ax_rho.set_ylabel(r"$\rho(r,z)/\rho_c(z)$") #\; [g \cdot cm^{-3}]$")
ax_rho.set_xlabel(r"$r/R_{\rm vir}$")
ax_rho.legend(handles=leg,  loc = "lower left")
fig_rho.savefig("Plots/density_profiles.pdf", bbox_inches='tight')

ax_Ts.set_xlim(xvec[0],xvec[-1])
#ax_Ts.set_ylim(1e1,1e5)
#ax_Ts.loglog([xvec[0],xvec[-1]],[Tcmb0*(1.+z),Tcmb0*(1.+z)],":",label="CMB")
ax_Ts.set_ylabel(r"$T_S\; [K]$")
ax_Ts.set_xlabel(r"$r/R_{\rm vir}$")
ax_Ts.legend(handles=leg, loc = "lower left")
fig_Ts.savefig("Plots/Ts.pdf", bbox_inches='tight')

ax_tau.set_xlim(aloverRv[0],aloverRv[-1])
#ax_tau.set_ylim(1e-3,1e1)
ax_tau.set_ylabel(r"$\tau$")
ax_tau.set_xlabel(r"$\alpha/R_{\rm vir}$")
ax_tau.legend(handles=leg, loc = "upper right")
fig_tau.savefig("Plots/optical_depth.pdf", bbox_inches='tight')
