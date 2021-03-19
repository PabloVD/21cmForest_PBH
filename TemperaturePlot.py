""""
Temperature plot code
Author: Pablo Villanueva Domingo
Started: October 2018
Last update: March 2021
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Source.functions import *


MassPBH_vals = [1., 10.,100.]
fpbh_vals = [1.e-2,1.e-3]

lins = ["-","--",":"]
cols = ["c","m","g"]
widths = [2.,1.,0.5]
colcdm = "b"


#--- TEMPERATURES ---#

fig_T, (ax_T) = plt.subplots(1,1,sharex=True)

Tleg = []

for MM, MassPBH in enumerate(MassPBH_vals):
    T_fpbh = []
    for ff, fpbh in enumerate(fpbh_vals):
        ThermalHisF = 'ThermalHistory_Mpbh_{:.1e}_fpbh_{:.1e}_Xi_30_Tmin_1.000e+04_Rfmp_15_chiUV_1.00e+30_Nalpha_4.00e+03.dat'.format(MassPBH, fpbh)
        file = np.loadtxt("21cmFiles/"+ThermalHisF,unpack=True)
        ax_T.plot(file[0],file[2],color=cols[ff],linestyle=lins[MM],lw=1.)#,linewidth=widths[MM],alpha=0.5)

#cdmfile = 'ThermalHistory_Mpbh_1.0e+00_fpbh_1.0e-08_Xi_30_Tmin_1.000e+04_Rfmp_15_chiUV_1.00e+30_Nalpha_4.00e+03.dat'
#file = np.loadtxt("21cmFiles/"+cdmfile,unpack=True)
#ax_T.plot(file[0],file[2],color=colcdm,alpha=0.1,lw=3.)
ax_T.plot(file[0],Tk_ad(file[0]),color=colcdm,alpha=1,lw=1.)
#ax_T.plot(file[0],Tk_ad(file[0]),color=colcdm)

Tleg.append(Line2D([0], [0], color=colcdm, linestyle="-", lw=4, label=r"$\Lambda$CDM") )
for ff, fpbh in enumerate(fpbh_vals):
    Tleg.append( Line2D([0], [0], color=cols[ff], linestyle="-", lw=4, label=r"$f_{PBH}=$"+scinot(fpbh)) )
for MM, MassPBH in enumerate(MassPBH_vals):
    Tleg.append( Line2D([0], [0], color="k", linestyle=lins[MM], lw=2, label="$M_{PBH}=$"+scinot(MassPBH)+" $M_{\odot}$" ) )
ax_T.set_ylabel("T [K]")
ax_T.set_xlabel("z")
ax_T.set_yscale("log")
#ax_T.set_xlim([file[0][-1],file[0][0]])
ax_T.set_xlim(10.,file[0][0])
ax_T.set_ylim(2.,100.)
ticks = range(10,29,2)
ax_T.set_xticks(ticks)
ax_T.legend(handles=Tleg)
fig_T.savefig("Plots/temperatures_pbh.pdf", bbox_inches='tight')
plt.show()
