import sympy
import random
import numpy as np
import mpmath as mp
import math

from .simulation import Simulation
from .point3d import Point3D
from .layer import Layer

DEBUG = False

# Used for debuging - only use this for a single photon simulation, since the output becomes unreadable very fast.
# Should probably implement a proper logger.     
def debug(string):
    if DEBUG:
        print(f"{datetime.datetime.now()} : {string}")

class Photon:
    """
    A class that represents a photon in a Monte Carlo simulation of photon propagation through a tissue.
    
    Attributes:
        simulation (Simulation): An instance of the Simulation class that contains information about the tissue structure and optical properties.
        position (Point3D): A Point3D object that represents the current position of the photon in 3D space.
        direction (Point3D): A Point3D object that represents the current direction of the photon in 3D space.
        current_layer (int): An layer object that represents the current layer of the tissue the photon is in.
        alive (bool): A boolean value that indicates if the photon is still alive (i.e. has not been absorbed or exited the tissue).
        path_history (list): A list of Point3D objects that represent the path of the photon through the tissue.
        MFP (float): The mean free path of the photon, which is the distance it travels before being scattered or absorbed.
    
    Methods:
        tick: Advances the photon's position and direction based on the simulation parameters and the current optical properties of the tissue.
        step_length: Determines the step length of the photon, which is the distance it will travel in a single time step.
        calculate_endpoint: Calculates the endpoint of the photon's step.
        check_boundary_crossing: Checks if the photon will cross a boundary between tissue layers during its step.
        reflect_refract: Calculates the new direction of the photon after it has been reflected or refracted at a tissue boundary.
        update_position: Updates the photon's position to a new point.
        check_inside: Checks if the photon is still inside the tissue.
        absorb_or_scatter: Determine if the photon is absorbed or scattered in the current tissue layer.
        scatter: A method to scatter the photon by modifying its direction according to the Henyey-Greenstein phase function.
        get2D_path_history: A method to get the 2D path of the photon in (r, z) format. This is primarily used for plotting. 
        update_position: Update the position of the photon and add it to the path history.
        calculate_endpoint: Calculate the endpoint of the photon's path, which is the point where the photon would end up if it continues to propagate without any interactions.
        refraction_direction: Calculate the refraction direction of the photon after interacting with a surface with a given normal and refractive indices.
        find_intersection: Find the intersection between the path of the photon and a surface.
        calculate_distance: Calculate the distance between two points.
        angle_of_incidence: Calculate the angle of incidence between the direction of the photon and the surface normal.
        fresnel_coefficients: Calculate the Fresnel coefficients for a given angle of incidence and refractive indices.
        recalculate_MFP: Recalculates the mean free path (MFP) of the photon after it refracts from one layer to the next.
    """
    
    def __init__(self, simulation, position, direction, current_layer):
        """
        Initializes a Photon object.
        
        Args:
            simulation (Simulation): An instance of the Simulation class that contains information about the tissue structure and optical properties.
            position (Point3D): A Point3D object that represents the starting position of the photon in 3D space.
            direction (Point3D): A Point3D object that represents the starting direction of the photon in 3D space.
            current_layer (Layer): An Layer object that represents the starting layer of the tissue the photon is in.
        """
        self.simulation = simulation
        self.position = position 
        self.direction = direction
        self.current_layer = current_layer
        self.alive = True
        self.path_history = [position]  
        self.MFP = 0


    def evaluate(self):
        """
        Advances the photon's position and direction based on the simulation parameters and the current optical properties of the tissue.
        """
        while self.alive:
            self.tick()
        return self

    def tick(self):
        """
        Advances the photon's position and direction based on the simulation parameters and the current optical properties of the tissue.
        """
        debug("--- New Tick ---")
        self.MFP = self.step_length()
        debug("MFP:")
        debug(self.MFP)

        new_endpoint = self.calculate_endpoint()
        
        crossing_layer, intersection = self.check_boundary_crossing()
        while intersection:             
            debug("---- New Loop ---")
            debug("Intersection found at:")
            debug(intersection)
            
            # Skip if we are at zero position
            if (intersection == self.position.as_array()).all():
                debug("Skipping zero position")
                break
            
            #Update the MFP
            distance = self.calculate_distance(self.position, Point3D(*intersection))
            self.MFP -= distance
            debug("MFP")
            debug(self.MFP)
            
            # Move photon to the intersection point
            self.update_position(Point3D(*intersection))
            debug("Current Direction:")
            debug(self.direction.as_array())
               
            vector = self.reflect_refract(intersection, crossing_layer)
            
            self.direction = Point3D(*vector).normalized()
            debug("Refraction/Reflection Vector:")
            debug(self.direction.as_array())
            
            crossing_layer, intersection = self.check_boundary_crossing()
            new_endpoint = self.calculate_endpoint()
        
        # Updates the photon position
        if self.check_inside() and self.alive:
            self.update_position(new_endpoint)
        elif len(self.path_history) > 2: 
            self.alive = False
        else:
            self.update_position(new_endpoint)
        
        # Check if the photon is still inside the simulation area
        if self.alive:
            # Absorb or scatter
            self.absorb_or_scatter()

        return self
        
    
    def step_length(self):
        """
        Determines the length of the step through the mean free path length before an absorption or scattering event.

        Returns:
            float: The mean free path length
        """
        MFP = (- np.log(1 - random.uniform(0 , 1))) / \
              ((self.current_layer.scattering_coefficient + self.current_layer.absorption_coefficient))
        
        return MFP
    
    def check_inside(self):
        """
        Check if the photon is inside the simulation area.

        Returns:
            bool: True if the photon is inside, False otherwise.
        """
        boundary_z = self.simulation.tissue_model.layers[-1].z_range[0]
        
        boundary_r = self.simulation.sim_size
        r = math.sqrt(self.position.x**2 + self.position.y**2)
        
        if (self.position.z <= 0) or (self.position.z > boundary_z) or (r > boundary_r):
            return False
        else:
            return True
    
    def check_boundary_crossing(self):
        """
        Checks if there is a boundary crossing with the plane specified in a Layer object.

        Returns:
            tuple: A tuple of two values:
            (layer, intersection), where layer is the next Layer object if the photon has crossed the boundary,
            False otherwise. intersection is the intersection point of the photon path with the boundary plane.
        """
        
        layers = self.simulation.tissue_model.layers
        current_index = layers.index(self.current_layer)
        
        endpoint = self.calculate_endpoint()
        
        if endpoint.z > self.current_layer.z_range[1]:
            # If the next layer is reached return the next layer
            debug("Boundary crossed next")
            
            # Find intersection for next_layer - plane of this layer
            intersection = self.find_intersection(layers[current_index])
            
            # Only return next layer if we aren't in the last layer
            if current_index + 1 < len(layers):
                return layers[current_index + 1], intersection
            else:
                self.alive = False
            
        elif endpoint.z < self.current_layer.z_range[0]:
            # If the previous layer is reached return the previous layer
            debug("Boundary crossed previous")
            debug(self.position.as_array())
            debug(endpoint.as_array())
            debug(self.direction.as_array())
            
            # Find intersection for previous layer - plane of previous layer
            intersection = self.find_intersection(layers[current_index - 1])
            
            return layers[current_index - 1], intersection
        
        return False, False

    
    def reflect_refract(self, intersection, layer):
        """
        Determines if the photon is reflected or refracted at a boundary. Reflection occurs if a random number is less than
        or equal to the reflectance, which is calculated using the Fresnel equations (R = 1/2 * (r_TM**2 + r_TE**2).

        Args:
            intersection (numpy.array): The intersection point of the photon path with the boundary plane.
            layer (Layer): The layer where refraction or reflection is occuring.

        Returns:
            numpy.array: The direction vector of the photon after reflection or refraction.
        """
        # Here we calculate the direction of the normal vector by just checking if the photon is travelling up or down
        if self.direction.z >= 0:
            n = np.array([0,0,1])
        else:
            n = np.array([0,0,-1])
            
        angle = self.angle_of_incidence(self.direction.as_array(), n)

        n_1 = self.current_layer.refractive_index
        n_2 = layer.refractive_index

        # We need to take a look at the case were we have normal incidence (theta = pi/2)
        # This should take care of the division by zero error in fresnel_coefficients but doesn't
        # one should take a closer look at this.
        if abs(angle - (mp.pi/2)) < 1e-6:
            R = ((n_1-n_2)/(n_1+n_2))**2
        else:
            coeff = self.fresnel_coefficients(angle, n_1, n_2)
            debug("Fresnel Coefficients:")
            debug(coeff)
            R = 1/2 * (coeff[0]**2 + coeff[1]**2)
        
        RND = random.uniform(0, 1)
        debug("RND:")
        debug(RND)
        debug("R:")
        debug(R)
        
        # Reflection
        if(RND <= R):
            # Incidence Vector
            I = self.direction.normalized().as_array()
            # Normal Vector
            reflection_vector = I - 2  * I.dot(n) * n
            debug("Reflection")
            return reflection_vector
        # Refract
        else:
            self.recalculate_MFP(layer)
            refraction_vector = self.refraction_direction(self.direction, n, n_1, n_2)
            self.current_layer = layer
            debug("Refraction")
            # Making sure that if would enter the last layer we destroy the photon
            
            
            return refraction_vector
    
    def absorb_or_scatter(self):
        """
        A method to decide whether the photon is absorbed or scattered. 
        This is determined based on a random number generated between 0 and 1 and
        the ratio of scattering coefficient to the sum of scattering coefficient and
        absorption coefficient of the current layer.
        """
        RND = random.uniform(0, 1)
        ratio = self.current_layer.absorption_coefficient / (self.current_layer.scattering_coefficient + self.current_layer.absorption_coefficient )
        
        # Absorption
        if RND <= ratio: 
            self.alive = False
        else: 
            self.scatter()
    
    

    def scatter(self):
        """
        A method to scatter the photon by modifying its direction according to the Henyey-Greenstein phase function.
        Source: https://en.wikipedia.org/wiki/Monte_Carlo_method_for_photon_transport
        """
        g = self.current_layer.anisotropy_factor
        # Direction cosine
        ds = self.direction.normalized()
        muxs = ds.x
        muys = ds.y
        muzs = ds.z
        
        chi = random.uniform(0,1)
        if g == 0:
            mutheta = 1 - 2 * chi
        else:
            mutheta = (1/(2 * g)) * ((1 + g**2 - ((1-g**2)/(1-g+2*g*chi))**2))
        
        fi = 2 * np.pi * chi
        
        costheta = mutheta
        sintheta = np.sqrt(1.0 - costheta**2) # sin(theta)
        sinfi = np.sin(fi)
        cosfi = np.cos(fi)
        if muzs == 1.0:
            muxd = sintheta * cosfi
            muyd = sintheta * sinfi
            muzd = costheta
        elif muzs == -1.0:
            muxd = sintheta * cosfi
            muyd = -sintheta * sinfi
            muzd = -costheta
        else:
            denom = np.sqrt(1.0 - muzs**2)
            muzcosfi = muzs * cosfi
            muxd = sintheta * (muxs * muzcosfi - muys * sinfi) / denom + muxs * costheta
            muyd = sintheta * (muys * muzcosfi + muxs * sinfi) / denom + muys * costheta
            muzd = -denom * sintheta * cosfi + muzs * costheta

        self.direction = Point3D(muxd, muyd, muzd)
        
    def get2D_path_history(self):
        """
        A method to get the 2D path of the photon in (r, z) format.
        This is primarily used for plotting. 

        Returns:
            A numpy array of points representing the 2D path of the photon.
        """
        points = []
        
        for point in self.path_history: 
            x = float(point.x)
            y = float(point.y)
            z = float(point.z)

            r = np.sqrt(x**2 + y**2)
            
            if x < 0:
                r *= -1
            
            points.append([r, z])
            
        return np.array(points)
    
    def update_position(self, new_pos):
        """
        Update the position of the photon and add it to the path history.

        Args:
            new_pos (Point3D): The new position of the photon.
        """
        self.position = new_pos
        self.path_history.append(self.position)
        
    def calculate_endpoint(self):
        """
        Calculate the endpoint of the photon's path, which is the point where the photon
        would end up if it continues to propagate without any interactions.

        Returns:
            Point3D: The endpoint of the photon's path.
        """
        normalized_direction = self.direction.normalized()
        return Point3D(*(self.position.as_array() + normalized_direction.scaled(self.MFP, self.MFP, self.MFP).as_array()))
     
    def refraction_direction(self, direction, surface_normal, n1, n2):
        """
        Calculate the refraction direction of the photon after interacting with a surface
        with a given normal and refractive indices.
        
        Source: https://stackoverflow.com/questions/29758545/how-to-find-refraction-vector-from-incoming-vector-and-surface-normal

        Args:
            direction (Point3D): The direction of the incident photon.
            surface_normal (np.array): The normal of the surface that the photon is interacting with.
            n1 (float): The refractive index of the medium that the photon is coming from.
            n2 (float): The refractive index of the medium that the photon is entering.

        Returns:
            np.array: The refraction direction of the photon after the interaction.
        """
        incident_dir = direction.as_array()
        
        incident_dir = incident_dir / np.linalg.norm(incident_dir)
        surface_normal = -surface_normal / np.linalg.norm(surface_normal)
        
        # 
        n = n1 / n2
        cosI = - np.dot(surface_normal, incident_dir)
        sinT2 = n * n * (1.0 - cosI * cosI)
        if sinT2 > 1.0:
            return None # TIR
        cosT = np.sqrt(1.0 - sinT2)
        return n * incident_dir + (n * cosI - cosT) * surface_normal

    def find_intersection(self, layer):
        """
        Find the intersection between the path of the photon and a surface.

        Args:
            layer (Layer): The surface that the photon is interacting with.

        Returns:
            list: The point of intersection between the photon's path and the surface.
        """
        debug("Finding intersection ...")
        endpoint = self.calculate_endpoint()
        segment = sympy.Line3D(self.position.as_array(), endpoint.as_array())
        debug("Segment:")
        debug(segment.points)
        debug("Plane:")
        debug(layer.plane)
        intersection = segment.intersection(layer.plane)
        if intersection:
            return intersection[0].evalf(6).coordinates
    
    def calculate_distance(self, p1, p2):
        """
        Calculate the distance between two points.

        Args:
            p1 (Point3D): The first point.
            p2 (Point3D): The second point.

        Returns:
            float: The distance between the two points.
        """
        distance = ((p2.x - p1.x)**2 + (p2.y - p1.y)**2 + (p2.z - p1.z)**2) **0.5
        return distance
    
    def angle_of_incidence(self, dir_ratio, normal_vector):
        """
        Calculate the angle of incidence between the direction of the photon and the surface normal.

        Args:
            dir_ratio (np.array): The direction of the photon.
            normal_vector (np.array): The normal of the surface that the photon is interacting with.

        Returns:
            float: The angle of incidence in radians.
        """
        # Normalize the direction ratio and normal vector
        dir_ratio = dir_ratio / np.linalg.norm(dir_ratio)

        # Calculate the dot product
        dot_product = np.dot(dir_ratio, normal_vector)

        # Calculate the angle of incidence
        angle = np.arccos(dot_product)

        return angle
    
    def fresnel_coefficients(self, angle_of_incidence, n_1, n_2):
        """
        Calculate the Fresnel coefficients for a given angle of incidence and refractive indices.

        Args:
            angle_of_incidence (float): The angle of incidence between the photon and the surface.
            n_1 (float): The refractive index of the medium that the photon is coming from.
            n_2 (float): The refractive index of the medium that the photon is entering.

        Returns:
            list: A list containing the Fresnel coefficients for reflected s-polarization, reflected p-polarization,
                  transmitted s-polarization, and transmitted p-polarization, respectively.
        """
        # Total reflection
        if (n_1 >= n_2) and (angle_of_incidence > math.asin(n_2/n_1)):
            rp = 1
            rs = 1
            ts = 0
            tp = 0
        else:
        # Fresnel:
            try:
                refraction_angle = math.asin(n_1 * math.sin(angle_of_incidence) / n_2)

                rp = -(math.sin(angle_of_incidence - refraction_angle)/(math.sin(angle_of_incidence + refraction_angle)))
                rs = -(math.tan(angle_of_incidence - refraction_angle)/(math.tan(angle_of_incidence + refraction_angle)))

                tp = -(2* math.cos(angle_of_incidence) * math.sin(refraction_angle))/(math.sin(angle_of_incidence + refraction_angle))
                ts = (2 * math.cos(angle_of_incidence) * math.sin(refraction_angle))/(math.sin(angle_of_incidence + refraction_angle) * math.cos(angle_of_incidence - refraction_angle))
            # In an edgecase it happens that the angle_of_incidence and the refraction_angle are 0.0
            # for the first tick when the photon is sitting on the (0,0,0) boundary.
            # Therefore transmission of the photon happens because we don't really want any interaction at this first layer
            except ZeroDivisionError:
                rp = 0
                rs = 0
                ts = 1
                tp = 1

        return [rs, rp, ts, tp]
    
    def recalculate_MFP(self, next_layer):
        """
        Recalculates the mean free path (MFP) of the photon after it refracts from one layer to the next.
        
        The MFP is determined by the sum of the scattering and absorption coefficients of the current layer, 
        and is scaled by the ratio of the sum of the scattering and absorption coefficients of the current and 
        next layers, in order to account for changes in tissue properties.
        
        Args:
            next_layer (Layer): The layer that the photon is entering after it refracts.
        """
        debug("Recalculating MFP after refraction:")
        
        prev_MFP_scale = (self.current_layer.scattering_coefficient + self.current_layer.absorption_coefficient)
        next_MFP_scale = (next_layer.scattering_coefficient + next_layer.absorption_coefficient)
        ratio = next_MFP_scale / prev_MFP_scale
        
        self.MFP *= ratio