""""
Plot the bounds on the fraction of PBHs as DM from the 21 cm forest absorbers
Author: Pablo Villanueva Domingo
Started: March 2021
Last update: March 2021
"""
import matplotlib.pyplot as plt
from scipy import interpolate, ndimage
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator, NullFormatter
from Source.functions import *

# Plot contours for the upper limits on the abundance
def PlotBound(z, MassPBH_vals, fpbh_vals, ax):

    # Initialize array with a random big number
    bounds1 = np.full((len(MassPBH_vals), len(fpbh_vals)),1.e3)
    bounds2 = np.full((len(MassPBH_vals), len(fpbh_vals)),1.e3)

    for im, Mpbh in enumerate(MassPBH_vals):
        file = np.loadtxt("Outputs/num_absorbers_z_{:.0f}_Mpbh_{:.1e}_from_fpbh_{:.1e}_to_{:.1e}.txt".format(z, Mpbh,fpbh_vals[0],fpbh_vals[-1]),unpack=True)
        fpbhs, cond1, cond2 = file[0], file[1], file[2]

        for jf, fpbh in enumerate(fpbh_vals):
            bounds1[im,jf] = cond1[jf]
            bounds2[im,jf] = cond2[jf]

    #plt.imshow(np.transpose(bounds1))
    #plt.show()
    #exit()

    #xx, yy = np.meshgrid(MassPBH_vals, fpbh_vals)
    #ax.scatter(xx,yy,c="k",alpha=0.7)

    # Interpolate
    bounds1_int = interpolate.interp2d(MassPBH_vals, fpbh_vals, np.transpose(bounds1),kind="linear")
    bounds2_int = interpolate.interp2d(MassPBH_vals, fpbh_vals, np.transpose(bounds2),kind="linear")

    logMpbh_vec = np.linspace(np.log10(MassPBH_vals[0]), np.log10(MassPBH_vals[-1]))
    logfpbh_vec = np.linspace(np.log10(fpbh_vals[0]), np.log10(fpbh_vals[-1]))

    # Filtering
    smooth1 = ndimage.filters.gaussian_filter(bounds1_int(10**logMpbh_vec, 10**logfpbh_vec), 1.)
    smooth2 = ndimage.filters.gaussian_filter(bounds2_int(10**logMpbh_vec, 10**logfpbh_vec), 1.)

    #ax.contourf(logMpbh_vec, logfpbh_vec, smooth1, levels=[0.,1.e3], alpha=0.5)
    #ax.contourf(logMpbh_vec, logfpbh_vec, smooth2, levels=[0.,1.e3], alpha=0.5)
    ax.contour(10.**logMpbh_vec, 10.**logfpbh_vec, smooth1, levels=[0.], alpha=1., linestyles="solid", colors="purple")
    ax.contour(10.**logMpbh_vec, 10.**logfpbh_vec, smooth2, levels=[0.], alpha=1., linestyles="dotted", colors="purple")


if __name__ == "__main__":

    # Redshift for the bounds
    z = 10.

    fig, ax = plt.subplots()

    MassPBH_vals = [0.1,0.3, 1.,10.,100., 1000.]
    fpbh_vals = [1.e-4,1.e-3,1.e-2]
    PlotBound(z, MassPBH_vals, fpbh_vals, ax)

    #fpbh_vals = [1.e-2,5.e-2,1.e-1,5.e-1,1.]
    fpbh_vals = [1.e-3,1.e-2,5.e-2]
    MassPBH_vals = [0.1, 0.3, 1.]
    PlotBound(z, MassPBH_vals, fpbh_vals, ax)

    ax.set_ylabel(r"$f_{\rm PBH}$")
    ax.set_xlabel(r"$M_{\rm PBH} \; [M_\odot$]")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.grid(True, linestyle=":", zorder=1.e-2, which='major')

    fig.savefig("Plots/bounds_z_"+str(z)+".pdf", bbox_inches='tight')
    plt.show()
