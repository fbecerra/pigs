from include import *

#**********************************************************************
#                         DATA PARAMETERS
#**********************************************************************
start_data           = 01
end_data             = 226
colors               = ['b', 'k', 'r']
path                 = '/Volumes/DARK-MATTER/enzo/'
simulations          = ['g1e10s4e10dm1e12']
total_sim            = len(simulations)

#**********************************************************************
#                         DISK PARAMETERS
#**********************************************************************
n_bins               = 100
specified_center     = [0.5, 0.5, 0.5]
percentage_of_stars  = 0.95
disk_height          = 1
disk_radius          = 30.0

#**********************************************************************
#                   PARAMETERS
#**********************************************************************
LOG_PLOT             = 1
DENSITY_PDF          = 0
TEMPERATURE_PDF      = 0
MASS_PDF             = 1
PARTICLES_MASS_PDF   = 0

lognorm_pdf = lambda p,x : numpy.exp(-(numpy.log(x)-p[0])**2/(2*p[1]**2))/(x*p[1]*numpy.sqrt(2*numpy.pi))
#lognorm_pdf = lambda p,x : numpy.exp(-(numpy.log((x-p[2])/p[1])**2/(2*p[1]**2)))/((x-p[2])*p[1]*numpy.sqrt(2*numpy.pi))
norm_pdf = lambda p,x : numpy.exp(-(x-p[0])**2./(2.*p[1]**2.))/(p[1]*numpy.sqrt(2.*numpy.pi))

for sim in range(total_sim):

   snapshots = np.array([25,226])
   labels    = [r't$\sim$100 Myr', r't$\sim$1 Gyr']

   for data in snapshots:
      N     = numpy.zeros(n_bins, float)
      d_m   = numpy.zeros(n_bins, float)
   
      dir   = path+simulations[sim]+"/DD%04i/" %data
      file  = simulations[sim]+"_%04i" %data
      pf    = load(dir+file)
      dd = pf.h.disk(specified_center, [0, 0, 1],disk_radius/pf["kpc"], disk_height/pf["kpc"])

      if DENSITY_PDF:
         value = numpy.log10(dd["Density"]*pc**3/Msun)
      if MASS_PDF:
         value = numpy.log10(dd["CellMassMsun"])
      if PARTICLES_MASS_PDF:
         stars = (dd["particle_type"] == PARTICLE_TYPE_STAR)
         value = numpy.log10(dd["ParticleMassMsun"][stars])
      value = value.reshape(-1)
   
      if DENSITY_PDF:
         dmin, dmax = dd.quantities["Extrema"]("Density", non_zero=True, lazy_reader=False)[0]
         dmin, dmax = dmin*pc**3/Msun, dmax*pc**3/Msun
      if MASS_PDF:
         dmin, dmax = dd.quantities["Extrema"]("CellMassMsun", non_zero=True, lazy_reader=False)[0]
      if PARTICLES_MASS_PDF:
         dmin, dmax = np.min(dd["ParticleMassMsun"][stars]), np.max(dd["ParticleMassMsun"][stars])
      if LOG_PLOT:
         dmin, dmax = numpy.log10(dmin), numpy.log10(dmax)
      delta = (dmax-dmin)/n_bins
      for i in range(0,n_bins):
         d_i    = delta*i + dmin
         d_f    = delta*(i+1) + dmin
         d_m[i] = delta*(i+1./2.) + dmin
         if LOG_PLOT and not PARTICLES_MASS_PDF:
            d_i    = 10**d_i
            d_f    = 10**d_f
            d_m[i] = 10**d_m[i]
         if DENSITY_PDF:
            N[i]=float(len(numpy.where((d_i<=dd["Density"]*pc**3/Msun) & (dd["Density"]*pc**3/Msun<d_f))[0]))/float(len(value))
         if MASS_PDF:
            N[i]=float(len(numpy.where((d_i<=dd["CellMassMsun"]) & (dd["CellMassMsun"]<d_f))[0]))/float(len(value))
         if PARTICLES_MASS_PDF:
            if LOG_PLOT:
               N[i]=float(len(numpy.where((d_i<=numpy.log10(dd["ParticleMassMsun"])) & (numpy.log10(dd["ParticleMassMsun"])<d_f))[0]))/float(len(value))
            else:
               N[i]=float(len(numpy.where((d_i<=dd["ParticleMassMsun"]) & (dd["ParticleMassMsun"]<d_f))[0]))/float(len(value))
   
      pl.step(np.log10(d_m),numpy.log10(N), color=colors[np.where(snapshots == data)[0]]) #, label=labels[np.where(snapshots == data)[0]])

      # Manual LN fit to density PDF
      if DENSITY_PDF and data == 226:
         pl.plot(numpy.arange(-8.0,4.0, 0.25), numpy.log10(norm_pdf([-2.5,1.0], numpy.arange(-8.0,4.0,0.25)))-1.2, 'r--')  #g1e10s4e10dm1e12

      # Power-law fit to high-mass PDF
      if MASS_PDF and data == 226:
         high_mass=numpy.where((numpy.log10(d_m) > 4.0) & (numpy.log10(d_m) < 6.5))[0]
         polycoeffs = sc.polyfit(numpy.log10(d_m[high_mass]),numpy.log10(N[high_mass]), 1)
         pl.plot(numpy.array([2.0,8.0]),polycoeffs[1]+polycoeffs[0]*numpy.array([2.0,8.0]), 'r--', label=r'$\alpha$ ='+"{0:.3}".format(polycoeffs[0]))
         pl.legend(loc=2, prop={'size':22})
   
      if PARTICLES_MASS_PDF:
         pl.ylabel(r'log (dN/dlogM)')
         pl.xlabel(r'log M [$M_\odot$]')
         high_mass=numpy.where((d_m > 4.89) & (d_m < 5))[0]
         polycoeffs = sc.polyfit(d_m[high_mass],numpy.log10(N[high_mass]), 1)
         pl.plot(numpy.array([4.89,5]),polycoeffs[1]+polycoeffs[0]*numpy.array([4.89,5]), 'r--', label=r'\alpha ='+"{0:.3}".format(np.abs(polycoeffs[0])))
         pl.legend(loc=2, prop={'size':32})
         pl.savefig('./profiles/particles_mass_pdf_'+simulations[sim]+'.png')

   if DENSITY_PDF:
      pl.ylim(-6.0,-1.0)
      pl.xlim(-8,2.5)
      pl.ylabel(r'log(N/N$_0$)')
      pl.xlabel(r'log $\rho_{gas}$ [M$_{\odot}$ pc$^{-3}$]')
      pl.savefig('./profiles/new_volume_density_pdf_vs_time_'+simulations[sim]+'.png')
   if MASS_PDF:
      pl.ylim(-6,-1)
      pl.xlim(-4,7)
      pl.ylabel(r'log(dN/dM)')
      pl.xlabel(r'log M [$M_\odot$]')
      pl.savefig('./profiles/new_mass_pdf_vs_time_'+simulations[sim]+'.png')
   pl.show()
   pl.close()
