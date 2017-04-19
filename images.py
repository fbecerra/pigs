from include import *

#**********************************************************************
#                         DATA PARAMETERS
#**********************************************************************
start_data           = 220
end_data             = 220
simulations          = ['g1e10s4e10dm1e12']
#simulations          = ['g1e10s4e10dm1e12', 'g1e10s4e10dm1e11', 'g1e10s1e10dm1e12', 'g1e10s1e11dm1e12']
path                 = '/Volumes/DARK-MATTER/enzo/'
total_sim            = len(simulations)

slice_lenght         = 30 # kpc
specified_center     = [0.5, 0.5, 0.5]
field                = "Density" #"Temperature"
zmax                 = 1.0e-20 
zmin                 = 1.0e-28
PARTICLE_TYPE_STAR   = 2
GRIDS                = 0
CIRCLE               = 0
PARTICLES            = 0
AGE_PLOT             = 0
RANDOM_DISKS         = 0
percentage_of_stars  = [0.90, 0.95, 1.00]

for sim in range(total_sim):

   for snap in range (start_data, end_data+1):
     
      dir   = path+simulations[sim]+"/DD%04i/" %(snap)
      file  = simulations[sim]+"_%04i" %(snap)
      pf = load(dir+file)

      # Projection
      pc = PlotCollection(pf)
      if snap == 200:
         use_colorbar=True
      else:
         use_colorbar=False
      pc.add_projection(field, 2, weight_field="Density", center=specified_center, use_colorbar=use_colorbar)
      if field == 'density':
        pc.set_cmap("gist_gray")
      elif field == 'temperature':
        pc.set_cmap("autumn") 
      pc.set_width(slice_lenght, 'kpc')
      pc.set_zlim(zmin, zmax)
      label = field.lower()
      if not RANDOM_DISKS and not PARTICLES:
        pc.save("./images/"+label+"_"+file, format="png", override=True)
      if RANDOM_DISKS:
         d = {'linestyle':'dashed', 'edgecolor':'yellow', 'fill':False}
         fin=open("./data/random_circles_01.dat", "r")
         for line in fin:
            line=line.split()
            x_random=float(line[1])
            y_random=float(line[2])
            r_random=float(line[3])
            specified_center_random=[x_random, y_random, 0.0]
            pc.plots[-1].modify['sphere'](specified_center_random, r_random, circle_args=d)
         fin.close()
         label += "_random_disks"
         if not PARTICLES:
           pc.save("./images/"+label+"_"+file, format="png", override=True)

#      # Case of centered disks
#      if CENTERED_RINGS:
#         dd    = pf.h.all_data()
#         min_radius_code, max_radius_code = dd.quantities["Extrema"]("ParticleDiskRadiusCode", non_zero=True, lazy_reader=False)[0]
#         min_radius, max_radius = dd.quantities["Extrema"]("ParticleDiskRadiuskpc", non_zero=True, lazy_reader=False)[0]
#         delta=(max_radius_code-min_radius_code)/total_points
#         for index in range(total_points):
#            r2=min_radius_code+(index+1)*delta
#            pc.plots[-1].modify['sphere'](specified_center, r2, circle_args=d)

      if PARTICLES:
         dd    = pf.h.all_data()
         stars = (dd["particle_type"] == PARTICLE_TYPE_STAR)
         star_radius = sorted(dd["ParticleDiskRadiusCode"][stars])
         number_of_particles = len(star_radius)
         if number_of_particles != 0:
            pc.set_zlim(1.0e-10, 1.0e-9)
            if AGE_PLOT:
              max_age=pf.current_time-dd["creation_time"][stars].min()
              age_fraction=(pf.current_time-dd["creation_time"][stars])/max_age
              color=[]
              for i in range(len(age_fraction)):
                 color.append((age_fraction[i],0.0,1-age_fraction[i]))
            else:
              color='r'
            pc.plots[0].modify["particles"](1.0, col=color, ptype=PARTICLE_TYPE_STAR)
            if AGE_PLOT:
              label += "_particle_age"
              pc.save("./images/"+label+"_"+file, format="png", override=True)
            else:
              label += "_particle"
              pc.save("./images/"+label+"_"+file, format="png", override=True)
            if CIRCLE:
               rpart_max = star_radius[int(math.floor(percentage_of_stars[0]*number_of_particles))]
               d = {'linestyle':'dashed', 'edgecolor':'blue', 'fill':False}
               pc.plots[-1].modify['sphere'](specified_center, rpart_max, circle_args=d)
               rpart_max = star_radius[int(math.floor(percentage_of_stars[1]*number_of_particles))]
               d = {'linestyle':'dashed', 'edgecolor':'green', 'fill':False}
               pc.plots[-1].modify['sphere'](specified_center, rpart_max, circle_args=d)
               rpart_min, rpart_max= dd.quantities["Extrema"]("ParticleDiskRadiusCode", non_zero=True, lazy_reader=False)[0]
               d = {'linestyle':'dashed', 'edgecolor':'yellow', 'fill':False}
               pc.plots[-1].modify['sphere'](specified_center, rpart_max, circle_args=d)
               label += "_circle"
               pc.save("./images/"+label+"_"+file, format="png", override=True)
      elif GRIDS:
         pc.plots[0].modify["grids"]()
         label += "_grids"
         pc.save("./images/"+label+"_"+file, format="png", override=True)
