import numpy as np
import matplotlib.pyplot as plt
import raytracer as rt
import scipy.optimize as op

def absolute(x, a, x0):
    return a*np.abs(x-x0)

def quadratic(x, a, x0, y0):
    return a*(x - x0)**2 + y0

def cubic(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

def hyperbola(x, a, b, h, k):
    return b*np.sqrt(1+(x-h)**2/a**2)+k

def propagation(ray):
    #Refraction
    lens.propagate_ray(ray)
        
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

  
    
#Paramaters (all length in mm)
z0=100
thickness=5 
K=0.02 #Magnitude of curveture
n1=1
n2=1.5168
ar=15
n=2 #Orientation: 1 is positve facing source, 2 is plane facing source
output_distance = 250

#Construct Lens and Screen
lens = rt.PlanoConvex(z0, thickness, K, n1, n2, ar, n)
screen = rt.OutputPlane(output_distance)
screenarray = rt.OutputPlaneArray(202, 210, 50)

#Add a single ray
# ray = rt.Ray(np.array([0,0,0]), np.array([0,0.1,1]))
# propagation(ray)

#Add a beam
plt.figure(dpi=300)
beam = rt.Beam(10)
for ray in beam.switch_on():        
    propagation(ray)


for surface in lens.surface():    
    z0=surface._z0 
    K=surface.curvature()
    ar=surface._ar
    print(K)
    
    #Inserting Lens onto yz plane
    if surface.curvature() == 0: #Plane
        zlens = np.array([z0, z0])
        ylens = np.array([0, ar])
    if surface.curvature() > 0:  #Positive 
        zlens = np.linspace(z0,-np.sqrt((1/K)**2-ar**2)+(z0+1/K), 200)
        ylens = np.zeros(len(zlens))
        ylens[1:] = np.sqrt((1/K)**2 - (zlens[1:]-(z0+1/K))**2)
        ylens[0] = 0
    if surface.curvature() < 0:  #Negative
        zlens = np.linspace(z0,+np.sqrt((1/K)**2-ar**2)+(z0+1/K), 200)
        ylens = np.zeros(len(zlens))
        ylens[1:] = np.sqrt((1/K)**2 - (zlens[1:]-(z0+1/K))**2)
        ylens[0] = 0
           
    #2D Plot on yz plane    
    plt.plot(zlens, ylens, color="C1", label='Lens surface')
    plt.plot(zlens, -ylens, color="C1")

plt.plot([output_distance, output_distance], [ar, -ar], label='Screen')

plt.xlabel("z /mm")
plt.ylabel("y /mm")
#plt.axis('scaled')
#plt.title("Ray Tracer on zy plane")
plt.grid()
plt.legend()
plt.show()

beam.initial_coordinates()
screen.output()
#screenarray.output()

#screenarray.z0()

#%%Spot Distribution
z, radius = beam.separation_distribution()
intercept = np.zeros(5)
err = np.zeros(5)

for i in range(5):
    x = z[3:-1]
    y = radius[3:-1,i]
    
    plt.plot(x, y, 'x', color='grey')
    
    guess_a = 0.1
    guess_x0 = float(x[(y==np.amin(y)).nonzero()])
    p0 = [guess_a, guess_x0]
    
    fit = op.curve_fit(absolute, x, y, p0=p0)
    x_fit = np.linspace(z[3], z[-2], 100)
    data_fit = absolute(x_fit, *fit[0])
   
    plt.plot(x_fit, data_fit)
    #plt.show()
    
    intercept[i] = fit[0][1]
    err[i] = np.sqrt(fit[1][1,1])
    
plt.title("Spot Radis Distribution (Plano-convex Singlet)")
plt.xlabel("z")
plt.ylabel("Radius")
#plt.ylim(0,1.5)
plt.show()

spread = 0.5*(np.amax(intercept) - np.amin(intercept))
print("Distribution")
print("f = %.2f ± %.2f" %(np.mean(intercept), spread))

#%%RMS Spot fit by quadratic and cubic 
# RMS Spot fit by approximation fit
# RMS_z, RMS_radius = beam.RMS_spot_radius()


# x = RMS_z[3:-1]
# y = RMS_radius[3:-1]
# plt.plot(x, y, 'x', color='grey')

# #Data fit
# x_fit = np.linspace(z[3], z[-2], 100)

# guess_a = 1
# guess_b = 1
# guess_c = 1
# guess_d = 1
# guess_x0 = float(x[(y==np.amin(y)).nonzero()])
# guess_y0 = np.amin(RMS_radius)

# quadratic_p0 = [guess_a, guess_x0, guess_y0] #2nd order approx
# cubic_p0 = [guess_a, guess_b, guess_c, guess_d] #3rd order approx

# quadratic_fit = op.curve_fit(quadratic, x, y, p0=quadratic_p0)
# quadratic_data_fit = quadratic(x_fit, *quadratic_fit[0])
# plt.plot(x_fit, quadratic_data_fit, label='1st Order')

# cubic_fit = op.curve_fit(cubic, x, y, p0=cubic_p0)
# cubic_data_fit = cubic(x_fit, *cubic_fit[0])
# plt.plot(x_fit, cubic_data_fit, label='2nd Order')

# plt.title("RMS Spot Radius (Plano-convex Singlet)")
# plt.xlabel("z")
# plt.ylabel("RMS Radius")
# plt.legend()
# plt.show()
###

#%%RMS Spot fit by hyperbola fit
RMS_z, RMS_radius = beam.RMS_spot_radius()

x = RMS_z[3:-1]
y = RMS_radius[3:-1]
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

plt.title("RMS Spot Radius (Plano-convex Singlet)")
plt.xlabel("z")
plt.ylabel("RMS Radius")
plt.show()

RMS_f = fit[0][2] #h
RMS_f_err = np.sqrt(fit[1][2,2])

f_radius = fit[0][3] + fit[0][1] #k + b
err_f_radius = np.sqrt(np.sqrt(fit[1][3,3])**2 + np.sqrt(fit[1][1,1])**2)

print("RMS spot radius:")
print("f = %.10f ± %.10f" %(RMS_f, RMS_f_err))
print("RMS spot radius at focus = %.10f ± %.10f" %(f_radius, err_f_radius))

