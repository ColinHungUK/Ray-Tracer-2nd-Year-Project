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
  
def Optimisation(p):    
    #Unpack parameters
    K1, K2 = p

    #Construct Lens and Screen
    lens = rt.BestFormLens(z0, thickness, K1, K2, n1, n2, ar)
    screen = rt.OutputPlane(output_distance)
    
    #Shoot beam
    beam = rt.Beam(10)
    for ray in beam.switch_on():        
        propagation(ray, lens, screen)
        
    RMS_radius = beam.RMS_spot_radius_at_screen()
    
    return RMS_radius

def propagation_with_plot(ray, lens, screen):
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
    plt.plot(z,y, "-", color='red')
  
    
#Fixed parameters (all length in mm) (global variables)
z0=100
thickness=5
n1=1
n2=1.5168
ar = 40
output_distance = 198

#Initial guess curvatures (|1/K| must be less than ar)
guess_K1 = 0.01
guess_K2 = -0.01
p0 = np.array([guess_K1, guess_K2])

#Boundary values of curvatures
K_min = -1/ar
K_max = 1/ar

#Minimising
fit = op.fmin_tnc(Optimisation, 
                  p0, 
                  approx_grad=True,
                  bounds=[(K_min, K_max), (K_min, K_max)])

#%%Optimised Result
K1 = fit[0][0]
K2 = fit[0][1]

print("Optimised Curvatures (mm^-1):", K1, K2)

#Plotting Result
lens = rt.BestFormLens(z0, thickness, K1, K2, n1, n2, ar)
screen = rt.OutputPlane(output_distance)

#Add a beam
beam = rt.Beam(10)
plt.figure(dpi=300)
for ray in beam.switch_on():        
    propagation_with_plot(ray, lens, screen)

edge = []
c = ['C1', 'C2']
for index, surface in enumerate(lens.surface()):    
    K = surface._curveture
    
    #Inserting Lens onto yz plane
    if surface.curvature() == 0: #Plane
        zlens = np.array([surface._z0, surface._z0])
        ylens = np.array([0, ar])
    if surface.curvature() > 0:  #Positive 
        zlens = np.linspace(surface._z0,-np.sqrt((1/K)**2-surface._ar**2)+(surface._z0+1/K), 200)
        ylens = np.zeros(len(zlens))
        ylens[1:] = np.sqrt((1/K)**2 - (zlens[1:]-(surface._z0+1/K))**2)
        ylens[0] = 0
    if surface.curvature() < 0:  #Negative
        zlens = np.linspace(surface._z0,+np.sqrt((1/K)**2-surface._ar**2)+(surface._z0+1/K), 200)
        ylens = np.zeros(len(zlens))
        ylens[1:] = np.sqrt((1/K)**2 - (zlens[1:]-(surface._z0+1/K))**2)
        ylens[0] = 0
        
    edge.append(zlens[-1])
           
    #2D Plot on yz plane    
    label = "K=" + str(round(K,4))
    plt.plot(zlens, ylens, color=c[index], label=label)
    plt.plot(zlens, -ylens, color=c[index])

plt.plot([output_distance, output_distance], [ar, -ar], label='Screen')

plt.xlabel("z /mm")
plt.ylabel("y /mm")
#plt.axis('scaled')
#plt.title("Ray Tracer on zy plane")
plt.grid()
plt.legend()
#plt.xlim(0,250)
plt.show()

beam.initial_coordinates()
screen.output()

print("At z =", output_distance, "mm, RMS spot radius =", beam.RMS_spot_radius_at_screen(), 'mm')

if edge[0] > edge[1]:
    print("Surface overlapped!")




