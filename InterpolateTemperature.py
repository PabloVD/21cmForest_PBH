""""
Script to interpolate the temperature from 21cmFAST simulations at a given redsfhit as a function of the PBH parameters
Author: Pablo Villanueva Domingo
Started: March 2021
Last update: March 2021
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import interpolate
import time
from Source.functions import *
from ImpactParam import MaxImpactParam
from scipy import interpolate, ndimage
import pickle
import matplotlib as mpl
mpl.rcParams.update({'font.size': 14})

# Isotherm which passes through f_pbh=1.
def fpbh_isotherm(Mpbh):
    a, b = -1.566, 0.
    return 10.**(b)*Mpbh**a

# Employ Rbf interpolation rather than standard interp2d
use_rbf = 0
# No astro heating: 1.e30
# Astro heating: 1.e55, 1.e56
zetax = 1.e30
#zetax = 1.e56

alphval = .8

MassPBH_vals = [1.0, 10., 100., 1000.]
#fpbh_vals = [1.e-5, 5.e-5, 1.e-4, 5.e-4, 1.e-3, 5.e-3, 1.e-2, 5.e-2, 1.e-1, 5.e-1, 1.]
fpbh_vals = [1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.]

zvec = [8., 10., 15.]
#zvec = [10.]

for z in zvec:

    tempmatrix = np.full((len(MassPBH_vals), len(fpbh_vals)),1.e3)
    # Cases including accretion
    for im, Mpbh in enumerate(MassPBH_vals):
        for jf, fpbh in enumerate(fpbh_vals):
            ThermalHisF = 'ThermalHistory_Mpbh_{:.1e}_fpbh_{:.1e}_Xi_30_Tmin_1.000e+04_Rfmp_15_chiUV_{:.2e}_Nalpha_4.00e+03.dat'.format(Mpbh, fpbh, zetax)
            file = np.loadtxt("21cmFiles/"+ThermalHisF,unpack=True)
            Temp = interpolate.interp1d(file[0], file[2], fill_value="extrapolate")
            tempmatrix[im,jf] = Temp(z)

    xx, yy = np.meshgrid(MassPBH_vals, fpbh_vals)
    #plt.scatter(xx,yy,c="k",alpha=0.7)

    # Interpolate
    if use_rbf:
        temp_int = interpolate.Rbf(np.log10(xx), np.log10(yy), np.log10(np.transpose(tempmatrix)))
    else:
        temp_int = interpolate.interp2d(np.log10(MassPBH_vals), np.log10(fpbh_vals), np.log10(np.transpose(tempmatrix)),kind="linear")

    # Save interpolator to file
    with open('Outputs/temperature_interpolator_z_{:.1f}_from_fpbh_{:.1e}_to_{:.1e}_zetax_{:.2e}.pkl'.format(z,fpbh_vals[0],fpbh_vals[-1],zetax), 'wb') as f:
        pickle.dump(temp_int, f)

    logMpbh_vec = np.linspace(np.log10(MassPBH_vals[0]), np.log10(MassPBH_vals[-1]),100)
    logfpbh_vec = np.linspace(np.log10(fpbh_vals[0]), np.log10(fpbh_vals[-1]),100)
    if use_rbf: logMpbh_vec, logfpbh_vec = np.meshgrid(logMpbh_vec, logfpbh_vec)

    #print( temp_int(np.log10(3.), np.log10(2.e-1)), temp_int(np.log10(5.), -1), temp_int(np.log10(20.), -2), temp_int(np.log10(90.), -3) )

    # Filtering
    smooth = ndimage.filters.gaussian_filter(temp_int(logMpbh_vec, logfpbh_vec), 1.)
    #smooth = temp_int(logMpbh_vec, logfpbh_vec)
    np.save("Outputs/temperature_interpolated_z_{:.1f}_zetax_{:.2e}".format(z,zetax),smooth)

    plt.contourf(10.**logMpbh_vec, 10.**logfpbh_vec, smooth, alpha=alphval)

    #plt.plot( [3., 90.], [2.e-1, 1.e-3] )
    plt.plot(10.**logMpbh_vec, fpbh_isotherm(10.**logMpbh_vec), "r--")

    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(10.**logfpbh_vec[0], 10.**logfpbh_vec[-1])
    plt.ylabel(r"$f_{\rm PBH} = \Omega_\mathrm{PBH}/\Omega_\mathrm{DM}$")
    plt.xlabel(r"$M_{\rm PBH} \; [M_\odot$]")
    plt.colorbar(label=r"log$_{10}(T_k/{\rm K})$")
    plt.title(r"$z=${:.0f}".format(z))
    plt.savefig("Plots/IGM_temperature_z_{:.1f}_zetax_{:.2e}.pdf".format(z,zetax), bbox_inches='tight')
    #plt.show()
    plt.close()
