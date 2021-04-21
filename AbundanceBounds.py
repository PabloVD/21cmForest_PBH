""""
Plot the bounds on the fraction of PBHs as DM from the 21 cm forest absorbers
Author: Pablo Villanueva Domingo
Started: March 2021
Last update: April 2021
"""
import matplotlib.pyplot as plt
from scipy import interpolate, ndimage
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator, NullFormatter
from matplotlib.lines import Line2D
from Source.functions import *
import matplotlib as mpl
mpl.rcParams.update({'font.size': 14})


# Plot contours for the upper limits on the abundance
def PlotBound(z, MassPBH_vals, fpbh_vals, zetax, ax, lw=1., alph=1., col="purple"):

    # Employ Rbf interpolation
    use_rbf = 1
    # Size of the smoothing filter
    sigma_filter = 1.5

    # Initialize array with a random big number
    bounds1 = np.full((len(MassPBH_vals), len(fpbh_vals)),1.e3)
    bounds2 = np.full((len(MassPBH_vals), len(fpbh_vals)),1.e3)

    for im, Mpbh in enumerate(MassPBH_vals):
        file = np.loadtxt("Outputs/num_absorbers_z_{:.0f}_Mpbh_{:.1e}_from_fpbh_{:.1e}_to_{:.1e}_zetax_{:.2e}.txt".format(z, Mpbh,fpbh_vals[0],fpbh_vals[-1],zetax),unpack=True)
        fpbhs, cond1, cond2 = file[0], file[1], file[2]

        for jf, fpbh in enumerate(fpbh_vals):
            bounds1[im,jf] = cond1[jf]
            bounds2[im,jf] = cond2[jf]

    xx, yy = np.meshgrid(MassPBH_vals, fpbh_vals)
    # This is for plotting the points of the grid
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

    # Heating from astrophysical sources
    # No astro heating: 1.e30
    # Astro heating: 1.e55, 1.e56
    zetaxs = [1.e30, 1.e55, 1.e56]
    cols = ["b", "purple", "r"]

    MassPBH_vals = np.logspace(-1, 3, 24)
    fpbh_vals = np.logspace(-5, 0 , 30)

    zvec = [10., 15.]
    lws = [1., 4.]
    alphs = [1., 0.5]

    for i, z in enumerate(zvec):

        fig, ax = plt.subplots()

        for j, zetax in enumerate(zetaxs):
            PlotBound(z, MassPBH_vals, fpbh_vals, zetax, ax, 1., 1., cols[j])
            #PlotBound(z, MassPBH_vals, fpbh_vals, zetax, ax, lws[i], alphs[i], cols[j])

        # Accretion line constraint
        #ax.plot(MassPBH_vals, 0.1*(MassPBH_vals/1.e3)**(-1.59), "c--")
        
        # Shot noise line constraint
        #ax.plot(MassPBH_vals, 1.e-2*MassPBH_vals**(-1.), "m--")

        ax.set_ylabel(r"$f_{\rm PBH}$")
        ax.set_xlabel(r"$M_{\rm PBH} \; [M_\odot$]")

        ax.set_xscale("log")
        ax.set_yscale("log")

        customleg, customlabs = [], []
        customleg.append(Line2D([0], [0], color="k", lw=2, linestyle="-", label=r"$0.01<\tau<0.03$"))
        customleg.append(Line2D([0], [0], color="k", lw=2, linestyle=":", label=r"$\tau>0.03$"))

        for j, zetax in enumerate(zetaxs):
            lab = r"$L_X/{\rm SFR}=$"+scinot(10.**(np.log10(zetax)-16))+r" erg s$^{-1}$M$_\odot^{-1}$yr"
            if zetax==1.e30:    lab = "No astrophysical heating"
            customleg.append(Line2D([0], [0], color=cols[j], lw=2, linestyle="-",label=lab))
        ax.legend(handles=customleg, loc="lower left",fontsize=10)

        ax.set_title(r"$z=${:.0f}".format(z))
        ax.set_ylim(1.e-5,1.)
        ax.grid(True, linestyle=":", zorder=1.e-2, which='major')

        fig.savefig("Plots/bounds_z_"+str(z)+".pdf", bbox_inches='tight')
        #fig.savefig("Plots/bounds_z_"+str(z)+"_zetax_{:.2e}.pdf".format(zetax), bbox_inches='tight')
    #plt.show()
