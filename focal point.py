import numpy as np
import matplotlib.pyplot as plt
import raytracer as rt
import scipy.optimize as op

def absolute(x, a, x0): #absolute linear function
    return a*np.abs(x-x0)

def hyperbola(x, a, b, h, k): #h, k are refered as i, j in summary report to avoid confusion
    return b*np.sqrt(1+(x-h)**2/a**2)+k

def propagation(ray):
    #Refraction
    n, intercept = lens.propagate_ray(ray)
    #normal1 =  n + intercept
    #normal2 = -n + intercept
    
    #Projection
    screenarray.propagate_ray(ray)
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
    #plt.plot([normal2[2], intercept[2], normal1[2]], [normal2[1], intercept[1], normal1[1]], "--", color='grey')
  
    
def propagation_check_focal(ray):
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
    
#%%    
#Paramaters (all length in mm)
z0 = 100
K=0.03
n1=1
n2=1.5
ar=1/0.03
output_distance = 250

#Construct Lens and Screen
lens = rt.SphericalRefraction(z0, K, n1, n2, ar)
screen = rt.OutputPlane(output_distance)
screenarray = rt.OutputPlaneArray(195, 203, 50)

#Add a single ray
# ray = rt.Ray(np.array([0,0,0]), np.array([0,0.1,1]))
# propagation(ray)

#Add a beam
beam = rt.Beam(10)
plt.figure(dpi=300.0)
for ray in beam.switch_on():        
    propagation(ray)


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

plt.xlabel("z /mm")
plt.ylabel("y /mm")
#plt.axis('scaled')
#plt.title("Ray Tracer on zy plane")
plt.grid()
plt.legend()
plt.show()

plt.figure(dpi=300.0)
beam.initial_coordinates()
plt.figure(dpi=300.0)
screen.output()
#screenarray.output()


#screenarray.z0()

#%%Spot Distribution
z, radius = beam.separation_distribution()
intercept = np.zeros(5)
err = np.zeros(5)

plt.figure(dpi=300.0)   

for i in range(5):
    x = z[3:-1]
    y = radius[3:-1,i]
    
    plt.plot(x, y, 'x', color='grey')
    
    guess_a = 0.1
    guess_x0 = float(x[(y==np.amin(y)).nonzero()])
    p0 = [guess_a, guess_x0]
    
    line = "N=" + str(i+1)
    fit = op.curve_fit(absolute, x, y, p0=p0)
    x_fit = np.linspace(z[3], z[-2], 100)
    data_fit = absolute(x_fit, *fit[0])
   
    plt.plot(x_fit, data_fit, label=line)
    #plt.show()
    
    intercept[i] = fit[0][1]
    err[i] = np.sqrt(fit[1][1,1])

#plt.title("Separation Distribution")
plt.xlabel("z /mm")
plt.ylabel("Separation /mm")
plt.ylim(0,1.2)
plt.legend()
plt.show()

spread = 0.5*(np.amax(intercept) - np.amin(intercept))
print("Distribution")
print("f = %.2f ± %.2f" %(np.mean(intercept), spread))


#%%RMS Spot fit by hyperbola fit
RMS_z, RMS_radius = beam.RMS_spot_radius()

x = RMS_z[3:-1]
y = RMS_radius[3:-1]
plt.figure(dpi=300.0)
plt.plot(x, y, 'x', color='grey')

#Data fit
guess_a = 1
guess_b = 0.0001
guess_h = float(x[(y==np.amin(y)).nonzero()])
guess_k = 0
p0 = [guess_a, guess_b, guess_h, guess_k]

fit = op.curve_fit(hyperbola, x, y, p0=p0)
x_fit = np.linspace(RMS_z[3], RMS_z[-2], 100)
data_fit = hyperbola(x_fit, *fit[0])
   
plt.plot(x_fit, data_fit)

#plt.title("RMS Spot Radius")
plt.xlabel("z /mm")
plt.ylabel("RMS Spot Radius /mm")
plt.ylim(0,0.4)
plt.show()

RMS_f = fit[0][2] #h
RMS_f_err = np.sqrt(fit[1][2,2])

f_radius = fit[0][3] + fit[0][1] #k + b
err_f_radius = np.sqrt(np.sqrt(fit[1][3,3])**2 + np.sqrt(fit[1][1,1])**2)

print("RMS spot radius:")
print("f = %.10f ± %.10f" %(RMS_f, RMS_f_err))
print("RMS spot radius at focus = %.10f ± %.10f" %(f_radius, err_f_radius))


#%%    
#Construct Lens and Screen
lens = rt.SphericalRefraction(z0, K, n1, n2, ar)
screen = rt.OutputPlane(RMS_f)

#Add a beam
beam = rt.Beam(10)
plt.figure(dpi=300.0)
for ray in beam.switch_on():        
    propagation_check_focal(ray)

plt.figure(dpi=300.0)
screen.output()

