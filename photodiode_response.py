import numpy as np

# constants
pixel_area = 2e-6 # in m²
pixel_width = np.sqrt(pixel_area) * 1000 # in mm
pixel_length = np.sqrt(pixel_area) * 1000 # in mm
spectral_file = 'csv/blue_spectra_poor.csv' # 10mA/cm2
pd_resp_file = 'csv/vemd5080x01.csv'
# pd_resp_file = 'csv/vemd5060x01.csv'

# photodiode geometry variables:
# Estimate area of PD to be 3x3mm, z distance is 20mm.
pd_width = 3
pd_length = 3
z = 20
# Centered between all 4 pixels (~10mm edge to edge) = 5mm
pd_center_loc = 5

# UDC fixture & Hamamatsu photodiode in contact with pixel:
# pd_resp_file = 'csv/S3588-08.csv'
# pd_width = 3
# pd_length = 30
# z = 0.1


# pull spectral data
arr = np.genfromtxt(spectral_file, delimiter=',', skip_header=1)
wavelength = arr[:,0] # in nanometers
radiance = arr[:,1] # [W/sr/m²]

# pull responsivity curve [A/W]
arr = np.genfromtxt(pd_resp_file, delimiter=',', skip_header=1)

# create responsivity curve with same wavelength range and step of spectral data
pd_resp = np.interp(wavelength, arr[:,0], arr[:,1], left=0)

# hemispherical lambertian emission profile
def integrand(thetas,phis):
    t, p = np.meshgrid(thetas, phis)
    f = np.cos(t)/2
    return f
thetas = np.deg2rad(np.arange(-90,91,1))
phis = np.deg2rad(np.arange(0,181,1))
s = integrand(thetas,phis)
# solid angle = pi for lambertian 
solid_angle = np.trapz(np.trapz(s, thetas, axis=1), phis, axis=0)

# photocurrent if all 2pi photons of the emission profile are captured
photo_current_buttcoupled = np.sum(solid_angle * radiance * pixel_area * pd_resp)
print(f"2pi photocurrent: {photo_current_buttcoupled:e}")

# pixel as pt source (middle) to near side of PD
x1 = pd_center_loc - (pixel_width/2) - (pd_length/2)
theta_1 = int(np.round(np.rad2deg(np.arctan(x1/z))))

# pixel as pt source (middle) to far side of PD
x2 = pd_center_loc - (pixel_width/2) + (pd_length/2)
theta_2 = int(np.round(np.rad2deg(np.arctan(x2/z))))

phi_1 = 0
y = pd_width/2
phi_2 = 2*int(np.round(np.rad2deg(np.arctan(y/z))))

thetas = np.deg2rad(np.arange(theta_1,theta_2+1))
phis = np.deg2rad(np.arange(phi_1,phi_2+1))
s = integrand(thetas,phis)
solid_angle = np.trapz(np.trapz(s, thetas, axis=1), phis, axis=0)
photo_current = np.sum(solid_angle * radiance * pixel_area * pd_resp)

print(f"simulated photocurrent: {photo_current:e}")
print(np.round(100*photo_current/photo_current_buttcoupled,3),'% of light enters PD')
