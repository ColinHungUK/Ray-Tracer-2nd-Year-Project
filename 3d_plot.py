import numpy as np
import matplotlib.pyplot as plt
import raytracer as rt

def propagation(ray):
    #Refraction
    n, intercept = lens.propagate_ray(ray)
    
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
        
    #3D Plotting  
    ax.plot3D(x,y,z, color='red')
    ax.scatter3D(x,y,z, color='red')


#Paramaters
z0 = 5
K=1/10
n1=1
n2=1.5
ar=7
output_distance = 50

#Construct Lens and Screen
lens = rt.SphericalRefraction(z0, K, n1, n2, ar)
screen = rt.OutputPlane(output_distance)

#3D Plotting
fig = plt.figure()
ax = plt.axes(projection="3d")    
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')   

#Add a beam
beam = rt.Beam(5)
for ray in beam.switch_on():        
    propagation(ray)
    
plt.show()

# %matplotlib inline 
# %matplotlib qt  
