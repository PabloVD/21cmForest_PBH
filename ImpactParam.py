""""
Script to compute the maximum impact parameter
Author: Pablo Villanueva Domingo
Started: October 2018
Last update: March 2021
"""

import matplotlib.pyplot as plt
import time
from Source.functions import *

tauvec = np.logspace(-3,0,num=numtau)
newtauvec = np.delete( np.delete(tauvec,0) ,-1)
dlogtau = np.log(tauvec[1])-np.log(tauvec[0])

lins = ["-","--",":"]

# Compute maximum impact parameter
def MaxImpactParam(z):

    logMvec = np.linspace(np.log(MJeans(z,Tk_ad(z))),np.log(Mmin(1.e4,z)),num=nummass)
    logimpparam = np.linspace(np.log(1.e-5),0.,num=numalf)
    Mvec, impparam = np.exp(logMvec), np.exp(logimpparam)

    maximpactmatrix = []
    for tau in tauvec:
        maximpactvec = []
        for M in Mvec:
            maximpact = 0.
            for impa in impparam:
                opdep = OptDepth21(M,z,con(M,z),impa*Rvir(M,z))
                if opdep>=tau:
                    maximpact = impa
            maximpactvec.append(Rvir(M,z)*maximpact)
        maximpactmatrix.append( maximpactvec)

    return np.array(maximpactmatrix)


if __name__ == "__main__":

    time_ini = time.time()

    fig_imp, (ax_imp) = plt.subplots(1,1)

    zvec = [10, 15]
    for iz, z in enumerate(zvec):

        imp = MaxImpactParam(z)
        np.save("Outputs/max_impact_param_z_{:.1f}".format(z)+"_nummass_"+str(nummass)+"_numtau_"+str(numtau)+"_numalf_"+str(numalf),imp)

        ax_imp.loglog( tauvec, imp[:,0], color="r", linestyle=lins[iz], label=r"$z=$"+str(z)+r", $M=$"+scinot(Mvec[0])+r" $M_{\odot}$" )
        ax_imp.loglog( tauvec, imp[:,-1], color="b", linestyle=lins[iz], label=r"$z=$"+str(z)+r", $M=$"+scinot(Mvec[-1])+r" $M_{\odot}$" )

    ax_imp.set_ylabel(r"$\alpha_{max}/R_{vir}$")
    ax_imp.set_xlabel(r"$\tau$")
    ax_imp.legend()
    fig_imp.savefig("Plots/impact_param_nummass_"+str(nummass)+"_numtau_"+str(numtau)+"_numalf_"+str(numalf)+".pdf", bbox_inches='tight')
    plt.close(fig_imp)

    print("Total minutes elapsed:",(time.time()-time_ini)/60.)
