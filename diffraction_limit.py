import numpy as np
import matplotlib.pyplot as plt
import raytracer as rt
import scipy.optimize as op

def propagation(ray, lens, screen):
    #Refraction
    lens.propagate_ray(ray)
        
    #Projection
    screen.propagate_ray(ray)
    
    #Generate points for plotting the ray
    x = []
    y = []
    z = []
    for i in range(len(ray.vertices())): 
        x.append(ray.vertices()[i][0])
        y.append(ray.vertices()[i][1])
        z.append(ray.vertices()[i][2])
        
    #Plotting on yz plane    
    # plt.plot(z,y, "-", color='red') 
  
def RMS_radius_vs_beam_diameter(beam_diameter):    
    #Construct Lens and Screen
    lens = rt.SphericalRefraction(z0, K, n1, n2, ar)
    screen = rt.OutputPlane(output_distance)
    
    #Shoot beam
    beam = rt.Beam(beam_diameter)
    for ray in beam.switch_on():        
        propagation(ray, lens, screen)
        
    RMS_radius = beam.RMS_spot_radius_at_screen()
    
    return RMS_radius  
    
#Fixed parameters (all length in mm) (global variables)
z0=100
K = 0.03
n1=1
n2=1.5
ar = 1/0.03
output_distance = 250

f = 198.39 #estimated by RMS spot radius and hyperbola fit 
wavelength = 588e-6 

list_of_beam_diameter = np.linspace(0,2,100)
RMS_radius = np.zeros(len(list_of_beam_diameter))
diffraction_scale = np.zeros(len(list_of_beam_diameter))

for i, D in enumerate(list_of_beam_diameter):
    RMS_radius[i] = RMS_radius_vs_beam_diameter(D)
    diffraction_scale[i] =  wavelength*f/D

#Plotting for comparison
plt.figure(dpi=300)
#plt.title("Diffraction Scale & RMS Spot Radius")    
plt.plot(list_of_beam_diameter, RMS_radius, label="RMS Spot Radius")
plt.plot(list_of_beam_diameter, diffraction_scale, label="Diffraction Scale")
plt.legend()
plt.xlabel("Beam Diameter /mm")
plt.ylabel("Limiting Factor")
plt.show()

#Calculating Airy Disk Radius from diffraction scale (= diffraction scale * 1.22)
# plt.title("Airy Disk Radius & RMS Spot Radius")
# plt.plot(list_of_beam_diameter, RMS_radius, label="RMS Spot Radius")
# plt.plot(list_of_beam_diameter, 1.22*diffraction_scale, label="Airy Disk Radius")
# plt.legend()
# plt.xlabel("Beam Diameter /mm")
# plt.ylabel("Limiting Radius /mm")
# plt.show()







