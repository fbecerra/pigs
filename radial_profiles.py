from include import *

#**********************************************************************
#                         DATA PARAMETERS
#**********************************************************************
start_data           = 01
path                 = '/Volumes/DARK-MATTER/enzo/'
#simulations          = ['g1e10s4e10dm1e12', 'g1e10s1e10dm1e10', 'g1e10s1e12dm1e12']
simulations          = ['g1e10s4e10dm1e12']
#snapshots            = [[226],[226],[226]]
snapshots            = [[25, 50, 100]]
colors               = ['b', 'k', 'r']
total_sim            = len(simulations)

#**********************************************************************
#                         DISK PARAMETERS
#**********************************************************************
disk_radius          = numpy.array(range(1,60))/3.
initial_radius       = disk_radius[0]
delta_radius         = initial_radius/2.
#final_radius         = disk_radius[len(disk_radius)-1]
final_radius         = 30.0
disk_height          = 5                  # kpc
specified_center     = [0.5, 0.5, 0.5]
percentage_of_stars  = 0.95

Q                    = numpy.zeros(len(disk_radius), float) 
sigma_c              = numpy.zeros(len(disk_radius), float)
V_rot                = numpy.zeros(len(disk_radius), float)
M_dyn                = numpy.zeros(len(disk_radius), float) 
Omega                = numpy.zeros(len(disk_radius), float)
M_rot                = numpy.zeros(len(disk_radius), float)
disk_sigma_gas       = numpy.zeros(len(disk_radius), float)
ring_vtan            = numpy.zeros(len(disk_radius), float)
ring_dispersion      = numpy.zeros(len(disk_radius), float)
ring_sigma_gas       = numpy.zeros(len(disk_radius), float)
enclosed_gas_mass    = numpy.zeros(len(disk_radius), float)

TO_PLOT              = 2   #254
ALL_TIMES            = 1
ALL_SIMS             = 0
PLOT_PROPERTY        = lambda prop: TO_PLOT >> prop & 1
TOOMRE               = 1
VROT                 = 2
OMEGA                = 3
MROT                 = 4
MASS                 = 5
SIGMA_GAS            = 6
Q_K89                = 7
ENC_MASS             = 8
ADD_THEO_MASSES      = 0

if ADD_THEO_MASSES:
   def gas_factor(R):
      Rd = 3.5           # kpc
      return math.exp(-R/Rd)*(math.exp(R/Rd)-R/Rd-1)
   
   def DM_factor(R, M):
      c = 12.0
      r_s = (M*Msun/(4.*PI*(200./3.)*c**3*1.8791e-29*0.67**2*(0.3*(1+0.1)**3+0.7)))**(1./3.)/kpc
      return (np.log(1+R/r_s)-(R/r_s)/(1+R/r_s))/(np.log(1+c)-c/(1+c))
   
   def stars_factor(a,b,R):
      return R**2*np.sqrt(R**2/(R**2+(a+np.sqrt(z**2+b**2))**2)**3+z**2*(a+np.sqrt(z**+b**2))**2/((z**2+b**2)*(R**2+(a+np.sqrt(z**+b**2))**2)**3))

if TO_PLOT != 0:

   for sim in range(total_sim):

     for snap in snapshots[sim]:
   
       #LOAD DATA
       dir   = path+simulations[sim]+"/DD%04i/" %snap
       file  = simulations[sim]+"_%04i" %snap
       pf    = load(dir+file)
       for i in range(len(disk_radius)):
          my_disk=pf.h.disk(specified_center, [0, 0, 1], disk_radius[i]/pf["kpc"], disk_height/pf["kpc"])
          if PLOT_PROPERTY(ENC_MASS):
            enclosed_gas_mass[i]=my_disk.quantities["TotalQuantity"]("CellMassMsun", lazy_reader=True)[0]                   # Msun
            if TO_PLOT == 2**ENC_MASS : continue
          ring=numpy.where((disk_radius[i] - delta_radius < my_disk["DiskRadiuskpc"]) & (my_disk["DiskRadiuskpc"] <= disk_radius[i] + delta_radius))
          V_rot[i]=numpy.mean(my_disk["DiskTangentialVelocity"][ring])                                             # cm/s
          if TO_PLOT == 2**VROT: continue

          Omega[i]=V_rot[i]/(disk_radius[i]*kpc)                                                                   # s
          if TO_PLOT == 2**OMEGA: continue
          M_dyn[i]=V_rot[i]**2*(disk_radius[i]*kpc)/Grav/Msun                                                      # Msun
          if TO_PLOT == 2**MASS: continue

          disk_area=na.pi*(disk_radius[i]*kpc)**2                                                                  # cm^2
          ring_area=na.pi*((disk_radius[i] + delta_radius)**2-(disk_radius[i] - delta_radius)**2)*kpc**2           # cm^2
          disk_gas_mass=my_disk.quantities["TotalQuantity"]("CellMassMsun", lazy_reader=True)[0]                   # Msun
          ring_gas_mass=my_disk["CellMassMsun"][ring]                                                              # Msun
          disk_sigma_gas[i]=disk_gas_mass*Msun/disk_area                                                           # g/cm^2
          if TO_PLOT == 2**SIGMA_GAS: continue

          M_rot[i]=na.pi**4*Grav**2.0*disk_sigma_gas[i]**3.0/4.0/Omega[i]**4/Msun                                  # Msun
          if TO_PLOT == 2**MROT: continue

          ring_sigma_gas[i]=sum(ring_gas_mass)*Msun/ring_area                                                         # g/cm^2
          ring_vtan[i]=sum(my_disk["DiskTangentialVelocity"][ring]*ring_gas_mass)/sum(ring_gas_mass)                  # cm/s
          ring_velocity_squared=sum(my_disk["DiskVelocitySquared"][ring]*ring_gas_mass)/sum(ring_gas_mass)            # cm/s
          ring_dispersion[i]=na.sqrt(ring_velocity_squared)-ring_vtan[i]                                              # cm/s
          if i > 0 :
             dVdR=(ring_vtan[i]-ring_vtan[i-1])/((disk_radius[i]-disk_radius[i-1])*kpc)
             kappa=na.sqrt(2*ring_vtan[i]**2/(disk_radius[i]*kpc)**2*(1+disk_radius[i]*kpc/ring_vtan[i]*dVdR))
             Q[i]=ring_dispersion[i]*kappa/(na.pi*Grav*ring_sigma_gas[i])
#             Q[i]=ring_dispersion[i]*Omega[i]/(na.pi*Grav*ring_sigma_gas[i])
#             sigma_c[i]=alpha*kappa*ring_dispersion[i]/3.36/Grav
          if TO_PLOT == 2**TOOMRE: continue

       
       #PLOT
       if PLOT_PROPERTY(Q_K89):
          dVdR0=(ring_vtan[1]-ring_vtan[0])/((disk_radius[1]-disk_radius[0])*kpc)
          kappa0=na.sqrt(2*ring_vtan[0]**2/(disk_radius[0]*kpc)**2*(1+disk_radius[0]*kpc/ring_vtan[0]*dVdR0))
          sigma_c[0]=alpha*kappa0*6e5/3.36/Grav
          star_radius=0
          for last_sim in range(snap-9, snap+1):
             dir   = path+simulations[sim]+"/DD%04i/" %(last_sim)
             file  = simulations[sim]+"_%04i" %(last_sim)
             pf    = load(dir+file)
             dd    = pf.h.all_data()
             stars = (dd["particle_type"] == PARTICLE_TYPE_STAR)
             s_radius = sorted(dd["ParticleDiskRadiuskpc"][stars])
             number_of_particles = len(s_radius)
             rpart_max = s_radius[int(math.floor(percentage_of_stars*number_of_particles))]
             star_radius+=rpart_max/10
          pl.xlabel(r'$R/R_{95}$')
          pl.ylim(ymax=1.0)
          pl.xlim(xmax=2.0)
          pl.ylabel(r'$\Sigma_g/\Sigma_c$')
          pl.plot(disk_radius/star_radius, ring_sigma_gas/sigma_c, label=simulations[sim], color=colors[sim])
          if not ALL_SIMS and not ALL_TIMES:
             pl.savefig('./radial_profiles/Sigma_crit_K89_'+simulations[sim]+'.png')
             pl.show()
             pl.close()
       if PLOT_PROPERTY(TOOMRE):
          dVdR0=(ring_vtan[1]-ring_vtan[0])/((disk_radius[1]-disk_radius[0])*kpc)
          kappa0=na.sqrt(2*ring_vtan[0]**2/(disk_radius[0]*kpc)**2*(1+disk_radius[0]*kpc/ring_vtan[0]*dVdR0))
          Q[0]=ring_dispersion[0]*kappa0/(na.pi*Grav*ring_sigma_gas[0])
          pl.xlabel(r'Radius [kpc]')
          pl.ylim(ymax=2)
          pl.ylabel(r'Q')
          pl.plot(disk_radius, Q, label=simulations[sim], color=colors[sim])
          if not ALL_SIMS and not ALL_TIMES:
             pl.savefig('./radial_profiles/Q_'+simulations[sim]+'.png')
             pl.show()
             pl.close()
       if PLOT_PROPERTY(VROT):
          pl.xlabel(r'Radius [kpc]', weight='bold')
          pl.ylim(ymax=300)
          pl.ylabel(r'$V_{rot}$ [km/s]')
          pl.plot(disk_radius, V_rot/1e5, label=simulations[sim], color=colors[sim])
          location=0
          if not ALL_SIMS and not ALL_TIMES:
             pl.savefig('./radial_profiles/Vrot_'+simulations[sim]+'.png')
             pl.show()
             pl.close()
       if PLOT_PROPERTY(OMEGA):
          pl.xlabel(r'Radius [kpc]')
          pl.ylabel(r'Omega [$s^{-1}]$')
          pl.plot(disk_radius, Omega, label=simulations[sim], color=colors[sim])
          if not ALL_SIMS and not ALL_TIMES:
             pl.savefig('./radial_profiles/Omega_'+simulations[sim]+'.png')
             pl.show()
             pl.close()
       if PLOT_PROPERTY(MROT):
          pl.xlabel(r'Radius [kpc]')
          #pl.ylim(ymax=3.0e7)
          pl.ylabel(r'$M_{rot}$ [$M_\odot$]')
          pl.plot(disk_radius, M_rot, label=simulations[sim], color=colors[sim])
          location=0
          if not ALL_SIMS and not ALL_TIMES:
             pl.savefig('./radial_profiles/Mrot_'+simulations[sim]+'.png')
             pl.show()
             pl.close()
       if PLOT_PROPERTY(ENC_MASS):
          pl.xlabel(r'Radius [kpc]')
          pl.ylabel(r'$M_{gas}$ [$M_\odot$]')
          pl.plot(disk_radius, enclosed_gas_mass, label=simulations[sim], color=colors[sim])
          location=0
          if not ALL_SIMS and not ALL_TIMES:
             pl.savefig('./radial_profiles/Mgas_'+simulations[sim]+'.png')
             pl.show()
             pl.close()
       if PLOT_PROPERTY(MASS):
          pl.xlabel(r'Radius [kpc]')
          pl.ylabel(r'Mass [$M_\odot$]')
          pl.plot(disk_radius, M_dyn, label='Dynamical Mass')
          if ADD_THEO_MASSES:
             gas_mass  =float(sim.split('_')[0].split('g')[1].split('s')[0])*gas_factor(numoy.array(disk_radius))
             star_mass =float(sim.split('_')[0].split('g')[1].split('s')[1].split('dm')[0])*stars_factor(a[sim], b[sim], numpy.array(disk_radius))
             dm_mass   =float(sim.split('_')[0].split('g')[1].split('s')[1].split('dm')[1])*DM_factor(numpy.array(disk_radius), float(sim.split('_')[0].split('g')[1].split('s')[1].split('dm')[1]))
             pl.plot(disk_radius, gas_mass, label='Gas Mass')
             pl.plot(disk_radius, star_mass, label='Star Mass')
             pl.plot(disk_radius, dm_mass, label='DM Mass')
             pl.legend(loc=2)
          if not ALL_SIMS and not ALL_TIMES:
             pl.savefig('./radial_profiles/Mass_'+simulations[sim]+'.png')
             pl.show()
             pl.close()
       if PLOT_PROPERTY(SIGMA_GAS):
          pl.xlabel(r'Radius [kpc]')
          pl.ylabel(r'$\Sigma_{gas}$ [g~cm$^{-2}]$')
          pl.plot(disk_radius, disk_sigma_gas, label=simulations[sim], color=colors[sim])
          if not ALL_SIMS and not ALL_TIMES:
             pl.savefig('./radial_profiles/SigmaGas_'+simulations[sim]+'.png')
             pl.show()
             pl.close()

     if ALL_TIMES and not ALL_SIMS:
       #pl.legend(loc=location)
       if PLOT_PROPERTY(Q_K89):
          pl.savefig('./radial_profiles/Sigma_crit_K89_time_'+simulations[sim]+'.png')
       if PLOT_PROPERTY(TOOMRE):
          pl.savefig('./radial_profiles/Q_time_'+simulations[sim]+'.png')
       if PLOT_PROPERTY(VROT):
          pl.savefig('./radial_profiles/Vrot_time_'+simulations[sim]+'.png')
       if PLOT_PROPERTY(OMEGA):
          pl.savefig('./radial_profiles/Omega_time_'+simulations[sim]+'.png')
       if PLOT_PROPERTY(MROT):
          pl.savefig('./radial_profiles/Mrot_time_'+simulations[sim]+'.png')
       if PLOT_PROPERTY(MASS):
          pl.savefig('./radial_profiles/Mass_time_'+simulations[sim]+'.png')
       if PLOT_PROPERTY(ENC_MASS):
          pl.savefig('./radial_profiles/Mgas_time_'+simulations[sim]+'.png')
       if PLOT_PROPERTY(SIGMA_GAS):
          pl.savefig('./radial_profiles/Sigma_gas_time_'+simulations[sim]+'.png')
       pl.show()
       pl.close()

if ALL_SIMS:
#   pl.legend(loc=location)
   if PLOT_PROPERTY(Q_K89):
      pl.savefig('./radial_profiles/Sigma_crit_K89_All.png')
   if PLOT_PROPERTY(TOOMRE):
      pl.savefig('./radial_profiles/Q_All.png')
   if PLOT_PROPERTY(VROT):
      pl.savefig('./radial_profiles/Vrot_All.png')
   if PLOT_PROPERTY(OMEGA):
      pl.savefig('./radial_profiles/Omega_All.png')
   if PLOT_PROPERTY(MROT):
      pl.savefig('./radial_profiles/Mrot_All.png')
   if PLOT_PROPERTY(MASS):
      pl.savefig('./radial_profiles/Mass_All.png')
   if PLOT_PROPERTY(ENC_MASS):
      pl.savefig('./radial_profiles/Mgas_All.png')
   if PLOT_PROPERTY(SIGMA_GAS):
      pl.savefig('./radial_profiles/Sigma_gas_All.png')
   pl.show()
   pl.close()
