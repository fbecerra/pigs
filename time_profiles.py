from include import *

#**********************************************************************
#                         DATA PARAMETERS
#**********************************************************************
start_data           = 1
every_data           = 15
end_data             = 226
colors               = ['b', 'k', 'r']
simulations          = ['g1e10s4e10dm1e12']
path                 = '/Volumes/DARK-MATTER/enzo/'
total_sim            = len(simulations)

#**********************************************************************
#                         DISK PARAMETERS
#**********************************************************************
disk_radius          = 30  #kpc
disk_height          = 5                  # kpc
specified_center     = [0.5, 0.5, 0.5]
total_points         = (end_data-start_data)/every_data + 1


TO_PLOT              = 4   #254
ALL_SIMS             = 1
PLOT_PROPERTY        = lambda prop: TO_PLOT >> prop & 1
SFR                  = 1
ENC_MASS             = 2
enclosed_gas_mass    = numpy.zeros(total_points, float) 
time                 = numpy.zeros(total_points, float)

for sim in range(total_sim):
  
   if PLOT_PROPERTY(SFR):
     dir   = path+simulations[sim]+"/DD%04i/" %end_data
     file  = simulations[sim]+"_%04i" %end_data
     pf    = load(dir+file)
     my_disk=pf.h.disk(specified_center, [0, 0, 1], final_radius/pf["kpc"], disk_height/pf["kpc"])
     sfr = StarFormationRate(pf, data_source=my_disk)
     pl.xlabel(r'Time [yr]')
     pl.ylabel(r'SFR [$M_\odot~yr^{-1}$]')
     pl.plot(sfr.time-sfr.time[0], sfr.Msol_yr, '-', linewidth=1.0, label=r'SFR$_{inst}$', color=colors[sim])
     pl.plot(sfr.time-sfr.time[0], sfr.Msol_cumulative/(sfr.time-sfr.time[0]), '--', label=r'SFR$_{avg}$', color=colors[sim])
     pl.legend(loc=1, prop={'size':20})
     pl.ylim(ymin=0)
     if not ALL_SIMS:
       location=1
       pl.savefig('./time_profiles/SFR_vs_time_'+simulations[sim]+'.png')
       pl.show()
       pl.close()
     if TO_PLOT == 2**SFR:
       continue
 
   for index in range(total_points):

     #LOAD DATA
     idx   = index * every_data + start_data
     dir   = path+simulations[sim]+"/DD%04i/" %idx
     file  = simulations[sim]+"_%04i" %idx
     pf    = load(dir+file)
     my_disk=pf.h.disk(specified_center, [0, 0, 1], disk_radius/pf["kpc"], disk_height/pf["kpc"])
     if PLOT_PROPERTY(ENC_MASS):
       enclosed_gas_mass[index] = my_disk.quantities["TotalQuantity"]("CellMassMsun", lazy_reader=True)[0]                   # Msun
       time[index] = pf.current_time
       if TO_PLOT == 2**ENC_MASS:
         continue
      
   time = (time - time[0]) / sec_per_year
   print time
   pl.xlabel(r'Time [yr]')
   if PLOT_PROPERTY(ENC_MASS):
     pl.plot(time, enclosed_gas_mass, color=colors[sim], label=simulations[sim])
     pl.ylabel(r'M$_{\rm enc}$ [M$_{\cdot}$]')
     if not ALL_SIMS:
       pl.savefig('./time_profiles/Enclosed_gas_mass_'+simulations[sim]+'.png')
       pl.show()
       pl.close()

if ALL_SIMS:
#   pl.legend(loc=location)
   if PLOT_PROPERTY(SFR):
      pl.savefig('./time_profiles/SFR_All.png')
   if PLOT_PROPERTY(ENC_MASS):
     pl.savefig('./time_profiles/Enclosed_gas_mass_All.png')
   pl.show()
   pl.close()
