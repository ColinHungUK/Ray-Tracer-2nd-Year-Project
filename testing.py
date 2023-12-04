import numpy as np
import matplotlib.pyplot as plt
import raytracer as rt

def propagation(ray): #Function to propagate ray
    #Refraction
    n, intercept = lens.propagate_ray(ray)
    normal1 =  n + intercept + 25
    normal2 = -n + intercept - 25
    
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
    plt.figure(dpi=300.0)
    plt.plot(z,y, "-", color='red') 
    plt.plot([normal2[2], intercept[2], normal1[2]], [normal2[1], intercept[1], normal1[1]], "--", color='grey')
    

#Paramaters
z0 = 100
K=-0.03
n1=1.5
n2=1
ar=1/0.03
output_distance = 150

#Construct Lens and Screen
lens = rt.SphericalRefraction(z0, K, n1, n2, ar)
screen = rt.OutputPlane(output_distance)

#Add a single ray
ray = rt.Ray(np.array([0,0,0]), np.array([0,0.1,1]))
propagation(ray)

#Add a beam
# beam = rt.Beam(10)
# for ray in beam.switch_on():        
#     propagation(ray)


#Inserting Lens onto yz plane
if K == 0: #Plane
    zlens = np.array([z0, z0])
    ylens = np.array([0, ar])
if K > 0:  #Positive 
    zlens = np.linspace(z0,-np.sqrt((1/K)**2-ar**2)+(z0+1/K), 200)
    ylens = np.zeros(len(zlens))
    ylens[1:] = np.sqrt((1/K)**2 - (zlens[1:]-(z0+1/K))**2)
    ylens[0] = 0
if K < 0:  #Negative
    zlens = np.linspace(z0,+np.sqrt((1/K)**2-ar**2)+(z0+1/K), 200)
    ylens = np.zeros(len(zlens))
    ylens[1:] = np.sqrt((1/K)**2 - (zlens[1:]-(z0+1/K))**2)
    ylens[0] = 0
    
    
#2D Plot on yz plane    
plt.plot(zlens, ylens, color="C1", label='Lens')
plt.plot(zlens, -ylens, color="C1")
plt.plot([output_distance, output_distance], [ar, -ar], label='Screen')


plt.xlabel("z")
plt.ylabel("y")
#plt.axis('scaled')
#plt.title("Ray Tracer on zy plane")
plt.grid()
plt.legend()
plt.show()

# beam.initial_coordinates()
#screen.output()


#%%3D Plot
# fig = plt.figure()
# ax = plt.axes(projection="3d")    
# ax.plot3D(x,y,z)
# ax.scatter3D(x,y,z)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')   
# plt.show()

