""""
Main script to compute the number of absorbers of the 21 cm forest
Author: Pablo Villanueva Domingo
Started: October 2018
Last update: March 2021
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import interpolate
import time
from Source.functions import *
from ImpactParam import MaxImpactParam


time_ini = time.time()


MassPBH_vals = [0.1, 0.3, 1., 10., 100., 1000.]
fpbh_vals = [1.e-3,1.e-2,5.e-2]
MassPBH_vals = [1., 10., 100., 1000.]
fpbh_vals = [1.e-3, 1.e-2]
zvec = [10, 15]

lins = ["-","--",":"]
cols = ["c","m","g"]
widths = [2.,1.,0.5]
colcdm = "b"

# 1 for using heating from 21cmFAST, 0 for assuming adiabatic cooling
use_21cmFAST = 1
# 1 for plotting
do_plots = 1
# 1 for plotting the derivative of the number of absorbers respect to tau
plot_derivative = 0

tauvec = np.logspace(-3,0,num=numtau)
newtauvec = np.delete( np.delete(tauvec,0) ,-1)
dlogtau = np.log(tauvec[1])-np.log(tauvec[0])



# Compute the cumulative number of absorbers
def CumulativeFunctions(z,fpbh,MassPBH,CDM,T,maximpactmatrix,logMvec,ax_cum=0,ax_der=0):
    cumulvec = []
    dercumulvec = []

    for tt, tau in enumerate(tauvec):
        maximpactvec = maximpactmatrix[tt]

        stepfactor = np.heaviside(logMvec - np.log(MJeans(z,T)),1.)
        if CDM:
            intcumul = integrate.simps( dndlnM_CDM(np.exp(logMvec),z)*np.pi*np.asarray(maximpactvec)**2.*stepfactor, logMvec )
            #intcumul = integrate.quad(lambda lnM: dndlnM_CDM(np.exp(lnM),z)*np.pi*np.asarray(maximpactvec)**2., logMvec[0],logMvec[-1])[0]
        else:
            intcumul = integrate.simps( dndlnM_PBH(np.exp(logMvec),z,fpbh,MassPBH)*np.pi*np.asarray(maximpactvec)**2.*stepfactor, logMvec )
            #intcumul = integrate.quad(lambda lnM: dndlnM_CDM(np.exp(lnM),z)*np.pi*np.asarray(maximpactvec)**2., logMvec[0],logMvec[-1])[0]
        cumulvec.append( (1.+z)**2.*drdz(z)*intcumul/1.e6/(MpcToCm/hlittle) )

    if do_plots:
        if CDM:
            col = colcdm
        else:
            col = cols[fpbh_vals.index(fpbh)]
        ax_cum.loglog( tauvec, cumulvec, linestyle=lins[zvec.index(z)], color=col )

        # Compute the derivative respect to tau
        if plot_derivative:
            for t, tau in enumerate(tauvec):
                if t!=0 and t!=len(tauvec)-1:
                    dercumulvec.append((cumulvec[t+1]-cumulvec[t-1])/2./dlogtau)
            ax_der.loglog( newtauvec, -np.array(dercumulvec), linestyle=lins[zvec.index(z)], color=col )

    Nabs_fit = interpolate.interp1d(tauvec, cumulvec)
    return Nabs_fit

# Main routine to compute the number of absorbers
def Absorbers_Routine(fpbh_vals, MassPBH, use_21cmFAST, do_plots, plot_derivative):

    if do_plots:
        if plot_derivative:
            fig_cum, (ax_cum, ax_der) = plt.subplots(2,1,sharex=True)
            fig_cum.subplots_adjust(hspace=0,wspace=0)
            fig_cum.subplots_adjust(bottom=0.1, right=0.5, top=0.9)
        else:
            fig_cum, ax_cum = plt.subplots(figsize=(5,5))
            ax_der=0
            #ax_cum.set_aspect('equal')
    else:
        ax_cum=0; ax_der=0

    if use_21cmFAST:
        sufix = "MassPBH_{:.1e}".format(MassPBH)
    else:
        sufix = "adiabatic"

    # Get temperatures
    #Tcdm = [Tk_ad(z) for z in zvec]
    T_fpbh = []

    if MassPBH<=1.:
        for ff, fpbh in enumerate(fpbh_vals):
            T_fpbh.append( [Tk_ad(z) for z in zvec] )

    else:
        # Get temperatures from 21cmFAST
        if use_21cmFAST:

            # CDM
            #cdmfile = 'ThermalHistory_Mpbh_1.0e+00_fpbh_1.0e-08_Xi_30_Tmin_1.000e+04_Rfmp_15_chiUV_1.00e+30_Nalpha_4.00e+03.dat'
            #file = np.loadtxt("21cmFiles/"+cdmfile,unpack=True)
            #Temp = interpolate.interp1d(file[0], file[2], fill_value="extrapolate")

            for ff, fpbh in enumerate(fpbh_vals):
                ThermalHisF = 'ThermalHistory_Mpbh_{:.1e}_fpbh_{:.1e}_Xi_30_Tmin_1.000e+04_Rfmp_15_chiUV_1.00e+30_Nalpha_4.00e+03.dat'.format(MassPBH, fpbh)
                file = np.loadtxt("21cmFiles/"+ThermalHisF,unpack=True)
                Temp = interpolate.interp1d(file[0], file[2], fill_value="extrapolate")
                T_fpbh.append( [Temp(z) for z in zvec] )
                #T_fpbh.append( [file[2][i] for i in [6,17]] )#[2,9,20]] ) # Take temperatures at redshifts 30,20,10

        else:
            for ff, fpbh in enumerate(fpbh_vals):
                T_fpbh.append( [Tk_ad(z) for z in zvec] )

    for iz, z in enumerate(zvec):

        # Import or compute the maximum impact parameter at redshift z
        impactparamfile = "Outputs/max_impact_param_z_{:.1f}".format(z)+"_nummass_"+str(nummass)+"_numtau_"+str(numtau)+"_numalf_"+str(numalf)+".npy"
        if os.path.exists(impactparamfile):
            maximpactmatrix = np.load(impactparamfile)
        else:
            print("Maximum impact parameter file absent. Computing...")
            maximpactmatrix = MaxImpactParam(z)
            np.save(impactparamfile, maximpactmatrix)

        # CDM
        logMvec = np.linspace(np.log(MJeans(z,Tk_ad(z))),np.log(Mmin(1.e4,z)),num=nummass)
        Nabs_cdm = CumulativeFunctions(z,0.,0.,1,Tk_ad(z),maximpactmatrix,logMvec,ax_cum,ax_der)
        Nabs_cdm_1 = Nabs_cdm(0.01) - Nabs_cdm(0.03)
        Nabs_cdm_2 = Nabs_cdm(0.03)
        print(z, "CDM", "{:.0f}".format(Nabs_cdm_1), "{:.0f}".format(Nabs_cdm_2))

        # PBH
        conds1, conds2 = [], []
        for ff, fpbh in enumerate(fpbh_vals):


            Nabs_pbh = CumulativeFunctions(z,fpbh,MassPBH,0,T_fpbh[ff][iz],maximpactmatrix,logMvec,ax_cum,ax_der)
            Nabs_pbh_1 = Nabs_pbh(0.01) - Nabs_pbh(0.03)
            Nabs_pbh_2 = Nabs_pbh(0.03)
            cond_1 = np.abs( Nabs_pbh_1-Nabs_cdm_1 ) - (np.sqrt(Nabs_pbh_1) + np.sqrt(Nabs_cdm_1))
            cond_2 = np.abs( Nabs_pbh_2-Nabs_cdm_2 ) - (np.sqrt(Nabs_pbh_2) + np.sqrt(Nabs_cdm_2))
            print(z, "{:.1e}".format(fpbh), "{:.1e}".format(MassPBH), "{:.0f}".format(Nabs_pbh_1), "{:.0f}".format(Nabs_pbh_2), "{:.0f}".format(cond_1), "{:.0f}".format(cond_2))
            conds1.append(cond_1); conds2.append(cond_2)

        np.savetxt("Outputs/num_absorbers_z_{:.0f}_Mpbh_{:.1e}_from_fpbh_{:.1e}_to_{:.1e}.txt".format(z, MassPBH,fpbh_vals[0],fpbh_vals[-1]),np.transpose([ fpbh_vals, conds1, conds2 ]),fmt=['%.1e','%.1f','%.1f'])

    if do_plots:
        customlegend = []
        customlegend.append( Line2D([0], [0], color=colcdm, lw=4, label=r"$\Lambda CDM$") )
        for ff, fpbh in enumerate(fpbh_vals):
            if use_21cmFAST:
                customlegend.append( Line2D([0], [0], color=cols[ff], linestyle="-", lw=4, label=r"$f_{PBH}=$"+scinot(fpbh)) )
            else:
                customlegend.append( Line2D([0], [0], color=cols[ff], linestyle="-", lw=4, label=r"$f_{PBH}M_{PBH}/M_{\odot}=$"+scinot(fpbh)) )
        for iz, z in enumerate(zvec):
            customlegend.append( Line2D([0], [0], color="k", linestyle=lins[iz], lw=2, label="$z=$"+str(z) ) )

        #ax_cum.set_ylim(1.e-2,1.e3)
        #ax_der.set_ylim(1.e-2,1.e3)
        ax_cum.set_ylabel(r"$dN(>\tau)/dz$")
        if plot_derivative:
            ax_der.set_ylabel(r"$|\tau d^2N/dzd\tau|$")
            ax_der.set_xlabel(r"$\tau$")
            ax_der.set_xlim(newtauvec[1],newtauvec[-1])
        else:
            ax_cum.set_xlabel(r"$\tau$")
            ax_cum.set_xlim(tauvec[0],tauvec[-1])
        ax_cum.legend(handles=customlegend, loc = "lower left", prop={"size":6})
        if use_21cmFAST:
            ax_cum.set_title(r"$M_{PBH}=$"+scinot(MassPBH)+"$ M_\odot$")


        fig_cum.savefig("Plots/absorbers_pbh_"+sufix+".pdf", bbox_inches='tight')
        #plt.show()
        plt.close(fig_cum)


#--- MAIN ---#

# Adiabatic case
Absorbers_Routine(fpbh_vals, 1., 0, do_plots, 1)

# Cases including accretion
for MassPBH in MassPBH_vals:
    Absorbers_Routine(fpbh_vals, MassPBH, use_21cmFAST, do_plots, plot_derivative)

print("Total minutes elapsed:",(time.time()-time_ini)/60.)
