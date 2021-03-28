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
def PlotBound(z, MassPBH_vals, fpbh_vals, zetax, ax, lw=1., alph=1., col="purple"):

    # Employ Rbf interpolation
    use_rbf = 1
    # Size of the smoothing filter
    sigma_filter = 3.

    # Initialize array with a random big number
    bounds1 = np.full((len(MassPBH_vals), len(fpbh_vals)),1.e3)
    bounds2 = np.full((len(MassPBH_vals), len(fpbh_vals)),1.e3)

    for im, Mpbh in enumerate(MassPBH_vals):
        file = np.loadtxt("Outputs/num_absorbers_z_{:.0f}_Mpbh_{:.1e}_from_fpbh_{:.1e}_to_{:.1e}_zetax_{:.2e}.txt".format(z, Mpbh,fpbh_vals[0],fpbh_vals[-1],zetax),unpack=True)
        fpbhs, cond1, cond2 = file[0], file[1], file[2]

        for jf, fpbh in enumerate(fpbh_vals):
            bounds1[im,jf] = cond1[jf]
            bounds2[im,jf] = cond2[jf]

    #plt.imshow(np.transpose(bounds1))
    #plt.show()
    #exit()

    xx, yy = np.meshgrid(MassPBH_vals, fpbh_vals)
    #ax.scatter(xx,yy,c="k",alpha=0.7)

    # Interpolate
    MassPBH_vals, fpbh_vals = np.array(MassPBH_vals), np.array(fpbh_vals)
    if use_rbf:
        bounds1_int = interpolate.Rbf(np.log10(xx), np.log10(yy), np.transpose(bounds1),kind="linear")
        bounds2_int = interpolate.Rbf(np.log10(xx), np.log10(yy), np.transpose(bounds2),kind="linear")
    else:
        bounds1_int = interpolate.interp2d(np.log10(MassPBH_vals), np.log10(fpbh_vals), np.transpose(bounds1),kind="linear")
        bounds2_int = interpolate.interp2d(np.log10(MassPBH_vals), np.log10(fpbh_vals), np.transpose(bounds2),kind="linear")

    logMpbh_vec = np.linspace(np.log10(MassPBH_vals[0]), np.log10(MassPBH_vals[-1]))
    logfpbh_vec = np.linspace(np.log10(fpbh_vals[0]), np.log10(fpbh_vals[-1]))
    if use_rbf: logMpbh_vec, logfpbh_vec = np.meshgrid(logMpbh_vec, logfpbh_vec)

    # Filtering
    smooth1 = ndimage.filters.gaussian_filter(bounds1_int(logMpbh_vec, logfpbh_vec), sigma_filter)
    smooth2 = ndimage.filters.gaussian_filter(bounds2_int(logMpbh_vec, logfpbh_vec), sigma_filter)

    ax.contour(10.**logMpbh_vec, 10.**logfpbh_vec, smooth1, levels=[0.], alpha=alph, linewidths=lw, linestyles="solid", colors=col)
    ax.contour(10.**logMpbh_vec, 10.**logfpbh_vec, smooth2, levels=[0.], alpha=alph, linewidths=lw, linestyles="dotted", colors=col)
    #plt.contour(10.**logMpbh_vec, 10.**logfpbh_vec, bounds1_int(logMpbh_vec, logfpbh_vec), levels=[0.], alpha=1., linestyles="solid", colors="purple")
    #plt.contour(10.**logMpbh_vec, 10.**logfpbh_vec, bounds2_int(logMpbh_vec, logfpbh_vec), levels=[0.], alpha=1., linestyles="dotted", colors="purple")



if __name__ == "__main__":

    # 1 for using heating from 21cmFAST, 0 for assuming adiabatic cooling
    use_21cmFAST = 1
    # 1 for plotting
    do_plots = 0
    # 1 for plotting the derivative of the number of absorbers respect to tau
    plot_derivative = 0
    # Heating from astrophysical sources
    # No astro heating: 1.e30
    # Astro heating: 1.e55, 1.e56
    zetax = 1.e30

    # NO ASTRO
    MassPBH_vals = np.logspace(-1, 3, 16)
    fpbh_vals = np.logspace(-5, 0 ,20)

    # ASTRO HEAT
    #MassPBH_vals = np.logspace(-1, 3, 10)
    #fpbh_vals = np.logspace(-5, 0 , 10)

    zvec = [10., 15.]
    #zvec = [15.]
    lws = [1., 4.]
    alphs = [1., 0.5]

    fig, ax = plt.subplots()

    for i, z in enumerate(zvec):
        PlotBound(z, MassPBH_vals, fpbh_vals, zetax, ax, lws[i], alphs[i], "b")

    zetax = 1.e55
    for i, z in enumerate(zvec):
        PlotBound(z, MassPBH_vals, fpbh_vals, zetax, ax, lws[i], alphs[i], "purple")

    zetax = 1.e56
    for i, z in enumerate(zvec):
        PlotBound(z, MassPBH_vals, fpbh_vals, zetax, ax, lws[i], alphs[i], "r")

    ax.set_ylabel(r"$f_{\rm PBH}$")
    ax.set_xlabel(r"$M_{\rm PBH} \; [M_\odot$]")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.grid(True, linestyle=":", zorder=1.e-2, which='major')

    fig.savefig("Plots/bounds_z_"+str(z)+"_zetax_{:.2e}.pdf".format(zetax), bbox_inches='tight')
    plt.show()
