from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"
from yt.mods import *
from yt.analysis_modules.star_analysis.api import *

#**********************************************************************
#                           NEW FIELDS
#**********************************************************************

def _DiskSurfaceDensity(field, data):
    return data["Density"]*data["dz"]
def _ConvertDiskSurfaceDensity(data):
    return data.convert("cm")
add_field("DiskSurfaceDensity", function=_DiskSurfaceDensity,
           convert_function=_ConvertDiskSurfaceDensity,  units=r"\rm{g}/\rm{cm^2}")

def _ParticleDiskRadius(field, data):
    center = data.get_field_parameter("center")
    DW = data.pf.domain_right_edge - data.pf.domain_left_edge
    radius = na.zeros(data["particle_position_x"].shape, dtype='float64')
    for i, ax in enumerate('xy'):
        r = na.abs(data["particle_position_%s" % ax] - center[i])
        radius += na.minimum(r, na.abs(DW[i]-r))**2.0
    na.sqrt(radius, radius)
    return radius
add_field("ParticleDiskRadiusCode", function=_ParticleDiskRadius,
          validators=[ValidateParameter("center")],
          particle_type=True,
          display_name = "Particle Disk Radius (code)")

def _DiskRadius(field, data):
    center = data.get_field_parameter("center")
    DW = data.pf.domain_right_edge - data.pf.domain_left_edge
    radius = na.zeros(data["x"].shape, dtype='float64')
    for i, ax in enumerate('xy'):
        r = na.abs(data[ax] - center[i])
        radius += na.minimum(r, na.abs(DW[i]-r))**2.0
    na.sqrt(radius, radius)
    return radius

def _ConvertDiskRadiusCGS(data):
    return data.convert("cm")
add_field("ParticleDiskRadius", function=_ParticleDiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiusCGS, units=r"\rm{cm}",
          particle_type = True,
          display_name = "Particle Disk Radius")
add_field("DiskRadius", function=_DiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiusCGS, units=r"\rm{cm}")

def _ConvertDiskRadiusMpc(data):
    return data.convert("mpc")
add_field("DiskRadiusMpc", function=_DiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiusMpc, units=r"\rm{Mpc}",
          display_name = "Disk Radius")
add_field("ParticleDiskRadiusMpc", function=_ParticleDiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiusMpc, units=r"\rm{Mpc}",
          particle_type=True,
          display_name = "Particle Disk Radius")

def _ConvertDiskRadiuskpc(data):
    return data.convert("kpc")
add_field("ParticleDiskRadiuskpc", function=_ParticleDiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiuskpc, units=r"\rm{kpc}",
          particle_type=True,
          display_name = "Particle Disk Radius")
add_field("DiskRadiuskpc", function=_DiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiuskpc, units=r"\rm{kpc}",
          display_name = "Disk Radius")

def _ConvertDiskRadiuskpch(data):
    return data.convert("kpch")
add_field("ParticleDiskRadiuskpch", function=_ParticleDiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiuskpc, units=r"\rm{kpc}/\rm{h}",
          particle_type=True,
          display_name = "Particle Disk Radius")
add_field("DiskRadiuskpch", function=_DiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiuskpc, units=r"\rm{kpc}/\rm{h}",
          display_name = "Disk Radius")

def _ConvertDiskRadiuspc(data):
    return data.convert("pc")
add_field("ParticleDiskRadiuspc", function=_ParticleDiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiuspc, units=r"\rm{pc}",
          particle_type=True,
          display_name = "Particle Disk Radius")
add_field("DiskRadiuspc", function=_DiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiuspc, units=r"\rm{pc}",
          display_name="Disk Radius")

def _ConvertDiskRadiusAU(data):
    return data.convert("au")
add_field("ParticleDiskRadiusAU", function=_ParticleDiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiusAU, units=r"\rm{AU}",
          particle_type=True,
          display_name = "Particle Disk Radius")
add_field("DiskRadiusAU", function=_DiskRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertDiskRadiusAU, units=r"\rm{AU}",
          display_name = "Disk Radius")
add_field("ParticleDiskRadiusCode", function=_ParticleDiskRadius,
          validators=[ValidateParameter("center")],
          particle_type=True,
          display_name = "Particle Disk Radius (code)")
add_field("DiskRadiusCode", function=_DiskRadius,
          validators=[ValidateParameter("center")],
          display_name = "Disk Radius (code)")

def obtain_rvec(data):
    center = data.get_field_parameter('center')
    coords = na.array([data['x'],data['y'],data['z']], dtype='float64')
    new_shape = tuple([3] + [1]*(len(coords.shape)-1))
    r_vec = coords - na.reshape(center,new_shape)
    return r_vec # axis 0 is the x,y,z

def obtain_velocities(data):
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = na.zeros(3, dtype='float64')
    xv = data["x-velocity"] - bv[0]
    yv = data["y-velocity"] - bv[1]
    zv = data["z-velocity"] - bv[2]
    return xv, yv, zv

def _DiskRadialVelocity(field, data):
    center = data.get_field_parameter("center")
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    new_field = ( (data['x']-center[0])*(data["x-velocity"]-bulk_velocity[0])
                + (data['y']-center[1])*(data["y-velocity"]-bulk_velocity[1])
                )/data["DiskRadiusCode"]
    if na.any(na.isnan(new_field)): # to fix center = point
        new_field[na.isnan(new_field)] = 0.0
    return new_field
def _DiskRadialVelocityABS(field, data):
    return na.abs(_DiskRadialVelocity(field, data))
def _ConvertDiskRadialVelocityKMS(data):
    return 1e-5
add_field("DiskRadialVelocity", function=_DiskRadialVelocity,
          units=r"\rm{cm}/\rm{s}",
          validators=[ValidateParameter("center")])
add_field("RadialVelocityABS", function=_DiskRadialVelocityABS,
          units=r"\rm{cm}/\rm{s}",
          validators=[ValidateParameter("center")])
add_field("DiskRadialVelocityKMS", function=_DiskRadialVelocity,
          convert_function=_ConvertDiskRadialVelocityKMS, units=r"\rm{km}/\rm{s}",
          validators=[ValidateParameter("center")])
add_field("DiskRadialVelocityKMSABS", function=_DiskRadialVelocityABS,
          convert_function=_ConvertDiskRadialVelocityKMS, units=r"\rm{km}/\rm{s}",
          validators=[ValidateParameter("center")])

def _ParticleDiskVelocityMagnitude(field, data):
    """M{|v|}"""
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    return ( (data["particle_velocity_x"]-bulk_velocity[0])**2.0 + \
             (data["particle_velocity_y"]-bulk_velocity[1])**2.0 )**(1.0/2.0)
add_field("ParticleDiskVelocityMagnitude", function=_ParticleDiskVelocityMagnitude,
          particle_type=True,
          take_log=False, units=r"\rm{cm}/\rm{s}")

def _DiskVelocityMagnitude(field, data):
    """M{|v|}"""
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    return ( (data["x-velocity"]-bulk_velocity[0])**2.0 + \
             (data["y-velocity"]-bulk_velocity[1])**2.0 )**(1.0/2.0)
add_field("DiskVelocityMagnitude", function=_DiskVelocityMagnitude,
          take_log=False, units=r"\rm{cm}/\rm{s}")

def _DiskVelocitySquared(field, data):
    return  data["DiskVelocityMagnitude"]**2
add_field("DiskVelocitySquared", function=_DiskVelocitySquared, units=r"\rm{cm}^2/\rm{s}^2")

def _DiskTangentialVelocity(field, data):
    return na.sqrt(data["DiskVelocityMagnitude"]**2.0
                 - data["DiskRadialVelocity"]**2.0)
add_field("DiskTangentialVelocity",
          function=_DiskTangentialVelocity,
          take_log=False, units=r"\rm{cm}/\rm{s}")

def _DiskAngularVelocity(field, data):
    r_vec = obtain_rvec(data)
    xv, yv, zv = obtain_velocities(data)
    diskr_vec = na.array([r_vec[0],r_vec[1]], dtype='float64')
    diskv_vec = na.array([xv,yv], dtype='float64')
    mom = na.cross(diskr_vec, diskv_vec, axis=0)
#    radsquare = r_vec[0]**2.0+r_vec[1]**2.0
#    momABS=na.sqrt(mom[0]**2.0+mom[1]**2.0+mom[2]**2.0)
    return mom/data["DiskRadiusCode"]**2
def _convertDiskAngularVelocity(data):
    return 1.0/data.convert("cm")
add_field("DiskAngularVelocity", function=_DiskAngularVelocity, convert_function=_convertDiskAngularVelocity,
           units=r"\rm{s}^{-1}", validators=[ValidateParameter('center')])
