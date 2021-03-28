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
import pickle
import matplotlib as mpl
from Source.functions import *
from ImpactParam import MaxImpactParam
from AbundanceBounds import PlotBound

mpl.rcParams.update({'font.size': 14})


# Some parameters
lins = ["-","--",":"]
cols = ["c","m","g"]
widths = [2.,1.,0.5]
colcdm = "b"

taumin1 = 0.01
taumin2 = 0.03

tauvec = np.logspace(-3,0,num=numtau)
newtauvec = np.delete( np.delete(tauvec,0) ,-1)
dlogtau = np.log(tauvec[1])-np.log(tauvec[0])


# Import the IGM temperature from interpolation or adiabatic
def GetTemperature(MassPBH, fpbh, z, use_21cmFAST, zetax, fpbh_vals, is_CDM=0):

    if use_21cmFAST:
        if is_CDM:
            MassPBH, fpbh = 1., 1.e-5

        # Take temperature from file if it exists
        ThermalHisF = 'ThermalHistory_Mpbh_{:.1e}_fpbh_{:.1e}_Xi_30_Tmin_1.000e+04_Rfmp_15_chiUV_{:.2e}_Nalpha_4.00e+03.dat'.format(MassPBH, fpbh, zetax)
        if os.path.exists("21cmFiles/"+ThermalHisF):
            file = np.loadtxt("21cmFiles/"+ThermalHisF,unpack=True)
            Temp = interpolate.interp1d(file[0], file[2], fill_value="extrapolate")
            Tk = Temp(z)
        # If file does not exist but it is above 1 Msun, interpolate
        elif MassPBH>=1.:
            #print("Using interpolation of temperature")
            with open('Outputs/temperature_interpolator_z_{:.1f}_from_fpbh_{:.1e}_to_{:.1e}_zetax_{:.2e}.pkl'.format(z,fpbh_vals[0],fpbh_vals[-1], zetax), 'rb') as f:
                temp_int = pickle.load(f)
            Tk = 10.**temp_int(np.log10(MassPBH), np.log10(fpbh))
        # Below 1 M_sun, temperature is like the CDM, neglect accretion
        else:
            Tk = GetTemperature(0., 0., z, use_21cmFAST, zetax, fpbh_vals, is_CDM=1)
    else:
        Tk = Tk_ad(z)

    return Tk


# Compute the cumulative number of absorbers
def CumulativeFunctions(z,fpbh,MassPBH,CDM,maximpactmatrix,logMvec,T_IGM,do_plots,plot_derivative,ax_cum=0,ax_der=0):
    cumulvec = []
    dercumulvec = []

    for tt, tau in enumerate(tauvec):
        maximpactvec = maximpactmatrix[tt]

        stepfactor = np.heaviside(logMvec - np.log(MJeans(z,T_IGM)),1.)
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
def Absorbers_Routine(zvec, fpbh_vals, MassPBH, use_21cmFAST, zetax, do_plots, plot_derivative):

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
        sufix = "MassPBH_{:.1e}_zetax_{:.2e}".format(MassPBH,zetax)
    else:
        sufix = "adiabatic"

    for iz, z in enumerate(zvec):

        # Import or compute the maximum impact parameter at redshift z
        impactparamfile = "Outputs/max_impact_param_z_{:.1f}".format(z)+"_nummass_"+str(nummass)+"_numtau_"+str(numtau)+"_numalf_"+str(numalf)+".npy"
        if os.path.exists(impactparamfile):
            maximpactmatrix = np.load(impactparamfile)
        else:
            print("Maximum impact parameter file absent. Computing...")
            maximpactmatrix = MaxImpactParam(z)
            np.save(impactparamfile, maximpactmatrix)

        logMvec = np.linspace(np.log(MJeans(z,Tk_ad(z))),np.log(Mmin(1.e4,z)),num=nummass)

        # CDM
        T_IGM = GetTemperature(0., 0., z, use_21cmFAST, zetax, fpbh_vals, 1)
        Nabs_cdm = CumulativeFunctions(z,0.,0.,1,maximpactmatrix,logMvec,T_IGM,do_plots,plot_derivative,ax_cum,ax_der)
        Nabs_cdm_1 = Nabs_cdm(taumin1) - Nabs_cdm(taumin2)
        Nabs_cdm_2 = Nabs_cdm(taumin2)
        print(z, "CDM", "{:.0f}".format(Nabs_cdm_1), "{:.0f}".format(Nabs_cdm_2))

        # PBH
        conds1, conds2 = [], []
        for ff, fpbh in enumerate(fpbh_vals):

            T_IGM = GetTemperature(MassPBH, fpbh, z, use_21cmFAST, zetax, fpbh_vals)
            Nabs_pbh = CumulativeFunctions(z,fpbh,MassPBH,0,maximpactmatrix,logMvec,T_IGM,do_plots,plot_derivative,ax_cum,ax_der)
            Nabs_pbh_1 = Nabs_pbh(taumin1) - Nabs_pbh(taumin2)
            Nabs_pbh_2 = Nabs_pbh(taumin2)
            cond_1 = np.abs( Nabs_pbh_1-Nabs_cdm_1 ) - (np.sqrt(Nabs_pbh_1) + np.sqrt(Nabs_cdm_1))
            cond_2 = np.abs( Nabs_pbh_2-Nabs_cdm_2 ) - (np.sqrt(Nabs_pbh_2) + np.sqrt(Nabs_cdm_2))
            print(z, "{:.1e}".format(fpbh), "{:.1e}".format(MassPBH), "{:.0f}".format(Nabs_pbh_1), "{:.0f}".format(Nabs_pbh_2), "{:.0f}".format(cond_1), "{:.0f}".format(cond_2))
            conds1.append(cond_1); conds2.append(cond_2)

        np.savetxt("Outputs/num_absorbers_z_{:.0f}_Mpbh_{:.1e}_from_fpbh_{:.1e}_to_{:.1e}_zetax_{:.2e}.txt".format(z, MassPBH,fpbh_vals[0],fpbh_vals[-1],zetax),np.transpose([ fpbh_vals, conds1, conds2 ]),fmt=['%.1e','%.1f','%.1f'])

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
        ax_cum.legend(handles=customlegend, loc = "lower left")#, prop={"size":6})
        if use_21cmFAST:
            ax_cum.set_title(r"$M_{PBH}=$"+scinot(MassPBH)+"$\, M_\odot$")


        fig_cum.savefig("Plots/absorbers_pbh_"+sufix+".pdf", bbox_inches='tight')
        #plt.show()
        plt.close(fig_cum)


#--- MAIN ---#
if __name__ == "__main__":

    time_ini = time.time()

    # 1 for using heating from 21cmFAST, 0 for assuming adiabatic cooling
    use_21cmFAST = 1
    # 1 for plotting
    do_plots = 0
    # 1 for plotting the derivative of the number of absorbers respect to tau
    plot_derivative = 0
    # 1 for computing bounds (equivalent to call AbundanceBounds.py separately)
    compute_bounds = 1
    # Heating from astrophysical sources
    # No astro heating: 1.e30
    # Astro heating: 1.e55, 1.e56
    zetax = 1.e30

    # FOR PLOTS
    MassPBH_vals = [1.,10.,100.,1000.]
    fpbh_vals = [1.e-3,1.e-2]

    # FOR COMPUTING BOUNDS (dense grid)
    MassPBH_vals = np.logspace(-1, 3, 16)
    fpbh_vals = np.logspace(-5, 0 ,20)

    # FOR COMPUTING BOUNDS (no so dense grid)
    #MassPBH_vals = np.logspace(-1, 3, 10)
    #fpbh_vals = np.logspace(-5, 0 , 10)

    zvec = [10., 15.]
    lws = [1., 4.]
    alphs = [1., 0.5]

    # Adiabatic case
    #Absorbers_Routine(fpbh_vals, 1., 0, do_plots, 1)

    # Cases including accretion
    for MassPBH in MassPBH_vals:
        print("Mass PBH: ",MassPBH,"\n")
        Absorbers_Routine(zvec, fpbh_vals, MassPBH, use_21cmFAST, zetax, do_plots, plot_derivative)

    # Compute the bounds on the PBH abundance
    if compute_bounds:

        fig, ax = plt.subplots()

        for i, z in enumerate(zvec):
            PlotBound(z, MassPBH_vals, fpbh_vals, zetax, ax, lws[i], alphs[i])

        ax.set_ylabel(r"$f_{\rm PBH}$")
        ax.set_xlabel(r"$M_{\rm PBH} \; [M_\odot$]")

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.grid(True, linestyle=":", zorder=1.e-2, which='major')

        fig.savefig("Plots/bounds_z_"+str(z)+"_zetax_{:.2e}.pdf".format(zetax), bbox_inches='tight')
        plt.show()

    print("Total minutes elapsed:",(time.time()-time_ini)/60.)
