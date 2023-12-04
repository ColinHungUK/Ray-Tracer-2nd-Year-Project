 # -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 09:47:45 2022

@author: Colin
"""
import numpy as np
import matplotlib.pyplot as plt

def normalise(r):
    '''
    Normalise a vector and return its unit vector
    
    r: the vector to be normalised
    
    mag: magnitude of the vector
    r/mag: unit vector
    
    '''
    mag = np.linalg.norm(r)
    return r/mag


def snells_law (i, n, n1, n2):
    '''
    Calculate refracted angle using Snell's Law
    
    i: incident ray direction vector
    t: transmitted/refracted ray direction vector
    n: surface normal, pointing toward n2, refracted side
    (ALL vectors are normalised)
    
    n1, n2: refractive index

    '''
    #Check for Total Internal Refraction
    sin_theta_i = np.linalg.norm(np.cross(i,n))/(np.linalg.norm(i)*np.linalg.norm(n))
    if sin_theta_i > n2/n1: 
        print("Total Interal Refraction")
        return None
    
    ratio = n1/n2 #Refractive index ratio
    
    # print("i, n, ratio", i, n, ratio)
    # print(1-(1-np.dot(n,i)**2)*ratio**2)
    
    parallel = np.sqrt(1-(1-np.dot(n,i)**2)*ratio**2)*n #component parallel to the normal vector
    perpendicular = ratio*(i - np.dot(n,i)*n) #component perpendicular to the normal vector
    # print("parallel=", parallel)
    # print("perpendicular=", perpendicular)
    
    t = parallel + perpendicular 
    
    #Snell's law check (good if both are equal)
    # sin_theta_t = np.linalg.norm(np.cross(t,n))/(np.linalg.norm(i)*np.linalg.norm(n))
    # print("n1/n2 =", ratio,"sin(theta_2)/sin(theta_1) =", sin_theta_t/sin_theta_i)
    
    return t



class Ray:
    """
    A class to represent an optical ray
    
    Positions and directions are represented by a 3D numpy array,
    stored together in a list
    
    k does not necessarily have to be normalised
    
    """
    def __init__(self, p=np.array([0,0,0]), k=np.array([0,0,0])):
        '''
        self._position: list of (1,3) numpy array,
                        store position vector of the ray
        self._direction: list of(1,3) numpy array,
                         store directional vector of the ray
        '''
        self._position = [p] 
        self._direction = [k]
        
    def p(self):
        '''
        Return the last position vector of the ray
        
        '''
        return self._position[-1]
    
    def k(self):
        '''
        Return the last position vector of the ray
        
        '''
        return self._direction[-1]
    
    def append(self, p, k):
        '''
        Append a new position vector and a direction vector to object
        
        '''
        self._position.append(p) #Position
        self._direction.append(k) #Direction        
        return self
    
    def vertices(self):
        '''
        Return all position vectors ray stored
        
        '''
        return self._position 
    

class Beam:
    '''
    An circular numpy array of parallel optical rays, representing a uniform collimated beam
    
    Contains objects of class 'Ray', but not a subclass of 'Ray'
    '''
    def __init__(self, R, nC=5):
        '''
        '''
        self._ray = [Ray([0,0,0],[0,0,1])]
        self._nC = nC
        r = R/nC
        self._nR = 6
        for i in range(nC):
            for j in range(self._nR*(i+1)):
                x = (r*(i+1)*
                     np.cos(2*np.pi/(self._nR*(i+1))*j))
                y = (r*(i+1)*
                     np.sin(2*np.pi/(self._nR*(i+1))*j))
                self._ray.append(Ray([x,y,0],[0,0,1]))               
                    
       
    def switch_on(self):
        '''
        Return every ray object stored in the beam object

        Returns
        -------
        self._ray: object of class Ray
            
        '''
        return self._ray
    
    def initial_coordinates(self):
        '''
        Plot the initial position of every Ray object in the xy plane

        Returns
        -------
        None.

        '''
        #Loop through every Ray object 
        for point in self._ray: 
            plt.plot(point.vertices()[0][0], point.vertices()[0][1], 'x', color="red")
        plt.xlabel("x /mm")
        plt.ylabel("y /mm")
        plt.axis('scaled')
        #plt.title("Beam")
        plt.grid()
        plt.show()
        
    def count(self):
        '''
        Give the number of Ray objects in the beam

        Returns
        -------
        len(self._ray): int
            Number of ray objects.

        '''
        return len(self._ray)
    
        
    def separation_distribution(self):
        '''
        Find the diameter of each ring of the beam at every position stored by ray objects

        Returns
        -------
        z : numpy array
            z position of every position vectors stored Ray.
        separation : numpy array
            separation between opposing points along the x-axis, for every ring of rays.

        '''
        n = len(self._ray[1].vertices())
        z = np.zeros(n)
        separation = np.zeros((n,self._nC))

        for i in range(n):
            for j in range(self._nC):
                a_index = int(self._nR/2*(j+1)**2 - self._nR/2*(j+1)) + 1                
                b_index = int(a_index + (self._nR*(j+1)/2))       
                #print(a_index, b_index)

                a = np.array(self._ray[a_index].vertices()[i])
                b = np.array(self._ray[b_index].vertices()[i])
                #print(a, b)
                separation[i,j] = np.linalg.norm(a-b)

            z[i] = self._ray[1].vertices()[i][2] #Return the z coordinates of one of the point
        
        return z, separation
    
    def RMS_spot_radius(self):
        '''
        Find the RMS spot radius at every position stored by the ray objects

        Returns
        -------
        z : numpy array
            z position of every position vectors stored Ray.
        RMSradius : numpy array
            RMS spot radius

        '''
        z, separation = self.separation_distribution()
        RMSradius = np.zeros(np.size(separation,0))
        
        #Deviation of each point from z-axis is given by
        #half the separation of the corresponding pair of points
        for i in range(np.size(separation,0)): #Loop through every z position
            sum_of_squared = 0
            count = 1 #Number of points considered, for mean calculation  
            for j in range(self._nC): #Loop through every rings of ray
                sum_of_squared += (separation[i,j]/2)**2*self._nR*(j+1)
                count += self._nR*(j+1)
            RMSradius[i] = np.sqrt(sum_of_squared/count)
                    
        return z, RMSradius
    
    def RMS_spot_radius_at_screen(self):
        '''
        Find the RMS spot radius at the last position ONLY

        Returns
        -------
        RMSradius : float
            RMS spot radius.

        '''
        z, radius = self.separation_distribution()
        RMSradius = np.zeros(np.size(radius,0))
        
        RMSradius = np.sqrt(np.mean(radius[-1,:]**2))
                    
        return RMSradius
        
  
class OpticalElement:
    '''
    Base class for optical elemtns (e.g. refracting surfaces, lenses)
    
    '''
    def propagate_ray(self, ray):
        raise NotImplementedError("Derived Class is missing the method")
        

class SphericalRefraction(OpticalElement):
    '''
    Class for spherical refracting surfaces.
    Derived from OpticalElement class.
    
    Method 'intercept' finds the position vector of the intercept.
    
    Assume that the z0 > 0
    and that ray is propagated form < z0
    
    '''
    def __init__(self, z0=0, K=0, n1=0, n2=0, ar=0):
        '''
        Initiating the class SphericalRefraction to form an lens object

        Parameters
        ----------
        z0 : float
            z-axis intersection of the surface toward the -z direction.
        thickness : float
            separation between the plane and the curved surface at the z-axis intersection.
        K : float
            Curveture of the positve surface.
        n1 : float
            Refractive index of the region outside of the lens.
        n2 : float
            Refractive index of the lens.
        ar : float
            Aperture radius.
        n : integer
            Set the orientation of the plano-convex.

        Raises
        ------
        Exception
            Aperture radius ar is greater than radius of curvature 1/K should the surface not be a plane.

        Returns
        -------
        None.
        '''
        self._z0 = z0 #intercept of the surface with the z-axis
        self._curveture = K #curvature of the surface
        self._n1 = n1 #refractive index on incident
        self._n2 = n2 #refractive index on refracted
        self._ar = ar #aperture radius
        if self._curveture != 0: #If the surface is not a plane
            self._O = np.array([0, 0, self._z0+1/self._curveture]) #Origin of the radius of curveture
            if ar > np.abs(1/K):
                raise Exception("Aperture radius is greater than radius of curvature")
        
    def intercept(self, ray):
        '''
        Find the intercept between a ray and the spherical refractive surface.

        Parameters
        ----------
        ray : Ray object
              
        Returns
        -------
        intercept: list
              Position vector of the intercept.

        '''
        unit_k = normalise(ray.k())
        
        #(1) Plane surface
        if self._curveture == 0:
            l = (self._z0 - ray.p()[2])/unit_k[2]
            intercept = [ray.p()[0]+l*unit_k[0],
                         ray.p()[1]+l*unit_k[1],
                         self._z0]
            return intercept
        
        #(2)Curved Surface
        r = ray.p() - self._O
        R = 1/self._curveture
        
        a = np.dot(r, r) - R**2
        b = np.dot(r, unit_k)
        
        #print("k=", ray.k())
        #print("unit_k=", unit_k)
        # print("b=", b, "a=", a)
        # print(self._O)
        
        if b**2 < a: #Complex roots -> no intercept
            print("No intercept")
            return None
        
        elif b**2 == a: #Repeated root -> 1 intercept
            l = -b**2 + np.sqrt(b**2 - a)  
            return l*unit_k + ray.p()
        
        elif b**2 > a: #Real roots -> 2 intercepts
            l1 = -b + np.sqrt(b**2 - a)
            l2 = -b - np.sqrt(b**2 - a)
            intercept1 = l1*unit_k + ray.p()
            intercept2 = l2*unit_k + ray.p()
                        
            #print(intercept1, intercept2)
        
            #Choose the intercept within the apperature radius
            if R > 0: #Positive
                offset = self._O[2]-np.sqrt(R**2 - self._ar**2)
                if offset >= intercept1[2] >= self._z0:
                    return intercept1
                elif np.isclose(intercept1[2], self._z0): #A cruel way to prevent floating point error
                    return intercept1
                elif np.isclose(intercept1[2], offset):
                    return intercept1 
                else:
                    return intercept2
                
            elif R < 0: #Negative
                offset = self._O[2]+np.sqrt(R**2 - self._ar**2)
                if offset <= intercept1[2] <= self._z0:
                    return intercept1
                elif np.isclose(intercept1[2], self._z0):
                    return intercept1
                elif np.isclose(intercept1[2], offset):
                    return intercept1 
                else:
                    return intercept2
            
    def propagate_ray(self, ray):
        '''
        Refract the ray based on the parameter of the surface using Snell's law. 

        Parameters
        ----------
        ray : Ray object
            The incident ray to be refracted.

        Returns
        -------
        n : (1,3) numpy array
            Normalised direction vector of the surface normal at the intercept.
        intercept : (1,3) numpy array
            Position vector of the intercept between the incident ray and the surface.

        '''
        intercept = self.intercept(ray)
        if self._curveture == 0: #Plane 
           n = normalise(np.array([0,0, intercept[2]])) 
        elif self._curveture > 0: #Positive
           n = normalise(self._O-intercept)
        elif self._curveture < 0: #Negative
           n = normalise(-self._O+intercept)
        
        #Find the refracted direction vector, 
        #by calling the Snell's law vector form function
        refracted = snells_law(ray.k(), 
                               n, 
                               self._n1, 
                               self._n2)
        
        #Append the new position and direction vectors at the intercept after refraction
        ray.append(intercept, refracted)
        
        return n, intercept
      
    def curvature(self):
        return self._curveture
    
    
class PlanoConvex(OpticalElement):
    '''
    Class representing a plano-convex singlet lens
    Consist of a convex refractive surface and a plane refractive surface
    
    '''
    def __init__(self, z0, thickness, K, n1, n2, ar, n):
        '''
        Initiating the class PlanoConvex to form an lens bject

        Parameters
        ----------
        z0 : float
            z-axis intersection of the surface toward the -z direction.
        thickness : float
            separation between the plane and the curved surface at the z-axis intersection.
        K : float
            Curveture of the positve surface.
        n1 : float
            Refractive index of the region outside of the lens.
        n2 : float
            Refractive index of the lens.
        ar : float
            Aperture radius.
        n : integer
            Set the orientation of the plano-convex.

        Raises
        ------
        Exception
            n must be either 1 or 2.

        Returns
        -------
        None.

        '''
        if n == 1: #positive surface facing the source
            self._element = [SphericalRefraction(z0=z0, K=K, n1=n1, n2=n2, ar=ar), 
                             SphericalRefraction(z0=z0+thickness, K=0, n1=n2, n2=n1, ar=ar)]
        elif n == 2: #plane facing the source
            self._element = [SphericalRefraction(z0=z0, K=0, n1=n1, n2=n2, ar=ar), 
                             SphericalRefraction(z0=z0+thickness, K=-K, n1=n2, n2=n1, ar=ar)]  
        else:
            raise Exception("n must be either 1 or 2")
    
    def propagate_ray(self, ray):
        '''
        Refract the ray based on the parameter of the surface using Snell's law. 

        Parameters
        ----------
        ray : Ray object
            The incident ray to be refracted.

        Returns
        -------
        n : (1,3) numpy array
            Normalised direction vector of the surface normal at the intercept.
        intercept : (1,3) numpy array
            Position vector of the intercept between the incident ray and the surface.

        '''
        for surface in self._element:
            n, intercept = surface.propagate_ray(ray)
        return n, intercept

    def surface(self):
        return self._element


class BestFormLens(OpticalElement):
    '''
    Class representing a dual surface lens.
    
    This class is used to optimising a best form lens,
    where curvatures of the two surfaces are optimised to minimise the RMS spot radius
    at a given focal point with a given apperature radius. 
    
    '''
    def __init__(self, z0, thickness, K1, K2, n1, n2, ar):
        self._element = [SphericalRefraction(z0=z0, K=K1, n1=n1, n2=n2, ar=ar), 
                         SphericalRefraction(z0=z0+thickness, K=K2, n1=n2, n2=n1, ar=ar)]
    
    def propagate_ray(self, ray):
        '''
        Refract the ray based on the parameter of the surface using Snell's law. 

        Parameters
        ----------
        ray : Ray object
            The incident ray to be refracted.

        Returns
        -------
        n : (1,3) numpy array
            Normalised direction vector of the surface normal at the intercept.
        intercept : (1,3) numpy array
            Position vector of the intercept between the incident ray and the surface.

        '''
        for surface in self._element:
            n, intercept = surface.propagate_ray(ray)
        return n, intercept

    def surface(self):
        return self._element
    
            
class OutputPlane(OpticalElement):
    '''
    Class for the output plane - a screen where the rays are projected at last.
    '''
    def __init__(self, z0):
        self._z0 = z0
        self._point = []
    
    def intercept(self, ray): 
        '''
        Find the intercept between a ray and the output plane

        Parameters
        ----------
        ray : Ray object
              
        Returns
        -------
        intercept: list
              Position vector of the intercept.

        '''
        unit_k = normalise(ray.k())
        l = (self._z0 - ray.p()[2])/unit_k[2]
        intercept = [ray.p()[0]+l*unit_k[0],
                     ray.p()[1]+l*unit_k[1],
                     self._z0]
        return intercept
    
    def propagate_ray(self, ray):
        '''
        Project the ray onto the output screen

        Parameters
        ----------
        ray : Ray object
            The incident ray.

        Returns
        -------
        n : (1,3) numpy array
            Normalised direction vector of the surface normal at the intercept.

        '''
        intercept = self.intercept(ray)
        n = np.array([intercept[0], intercept[1], self._z0]) 
        self._point.append(n)
        ray.append(intercept, ray.k())
        return n
    
    def output(self):
        '''
        Plot the positions of every ray projected onto the screen. 

        Returns
        -------
        None.

        '''
        for point in self._point:
            plt.plot(point[0], point[1], 'x', color='red')
        plt.xlabel("x /mm")
        plt.ylabel("y /mm")
        plt.axis('scaled')
        #title = "Output Plane at z=" + str(self._point[0][2])
        #plt.title(title)
        plt.grid()
        plt.show()
        return self._point
    
    def count(self):
        '''
        Give the number of ray projected onto the output plane

        Returns
        -------
        len(self._point): int
            Number of points on the output plane. 

        '''
        return len(self._point)
    

class OutputPlaneArray(OpticalElement):
    '''
    A composist class contain objects of OutputPlane.
    
    Represents an array of parallel transparent screens that do not interact with the ray.
    '''
    def __init__(self, zmin, zmax, n):
        interval = (zmax-zmin)/(n-1) #distance interval between screens, in mm
        self._plane = []
        for i in range(n):
            self._plane.append(OutputPlane(zmin+interval*i))
            
    def z0(self):
        for plane in self._plane:
            print(plane._z0)
          
    def propagate_ray(self, ray):
        '''
        Propagate the ray through the array of screens.
        Projected positions are appended to the ray objects and the plane

        Returns
        -------
        None.
        
        '''
        for plane in self._plane:
            intercept = plane.intercept(ray)
            n = np.array([intercept[0], intercept[1], plane._z0]) 
            plane._point.append(n)
            ray.append(intercept, ray.k())
            
    def output(self):
        '''
        Plot the positions of the rays projected onto every plane

        Returns
        -------
        None.

        '''
        for plane in self._plane:
            for point in plane._point:
                plt.plot(point[0], point[1], 'x', color='red')
            plt.xlabel("x")
            plt.ylabel("y")
            plt.axis('scaled')
            title = "Output Plane at z=" + str(plane._point[0][2])
            plt.title(title)
            plt.grid()
            plt.show()
            
    
    
            
    
            
    
    
        
        