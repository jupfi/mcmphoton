import numpy as np

from .photon import Photon
from .point3d import Point3D

class Source:
    """
    The Source class represents a light source in the Monte Carlo simulation of photon propagation in tissue.

    Attributes:
        simulation (Simulation): A reference to the parent Simulation instance.
        position (Point3D): The position of the source in 3D space.
        direction (Point3D): The direction in which the photons are emitted from the source.
        number_photons (int): The number of photons to be emitted from the source.
        photons (list): A list of Photon instances emitted from the source.
        photons_alive (list): A list of Photon instances that are still propagating in the tissue.

    Methods:
        create_photons: Creates number_photons Photon instances with the same direction as the source.
        move_photons_to_boundary: Moves all photons to the boundary between the first and second layers.
        tick: Updates the position and state of all Photon instances that are still propagating in the tissue.
        get2D_position: Converts the 3D position of the source to a projection in the (r, z) plane.
    """
    
    def __init__(self, simulation, position, direction, number_photons):
        """
        Args:
            simulation (Simulation): Simulation in which the source is used.
            position (Point3D): Position of the source in 3D space.
            direction (Point3D): Emission direction of the source.
            number_photons (int): Number of photons created from the source.
        """
        self.simulation = simulation
        self.position = position
        self.direction = direction
        self.number_photons = number_photons
        
        self.photons = []
        
    def create_photons(self):
        """
        Creates photons from the source.
        """
        for i in range(self.number_photons): 
            #Gets the 'zero' layer
            current_layer = self.simulation.tissue_model.layers[0]
            
            # Creates the photons with start direction the same as for the source. 
            # Here one could implement more complicated emission patterns
            photon = Photon(self.simulation, self.position, self.direction, current_layer)
            self.photons.append(photon)
    
    # Right now this implementation only moves the photons to the zero position
    def move_photons_to_boundary(self):
        """
        Moves all photons to the boundary of the tissue model.
        """
        for photon in self.photons:
            photon.update_position(Point3D(0, 0, 0))
            photon.current_layer = self.simulation.tissue_model.layers[1]
    
    # Evaluates the photons in the simulation with multiple threads using multiprocessing
    def evaluate_photons(self):
        """
        Evaluates the photons in the simulation with multiple threads using multiprocessing.
        """
        from multiprocessing import Pool, cpu_count
        # Create a pool of processes
        pool = Pool(processes=cpu_count())
        
        # Evaluate the photons
        result = pool.map(Photon.evaluate, self.photons)
        
        # Close the pool
        pool.close()
        pool.join()
        self.photons = result

    # Tick method for every tick of the Monte Carlo simulation
    def tick(self):
        """
        Updates the position of the photons in each time step of the simulation.
        """

        for photon in self.photons:
            if photon.alive:
                photon.tick()
    
    # Converts the 3D position of the source to a projection in the (r, z) plane
    def get2D_position(self):
        """
        Converts the 3D position of the source to a projection in the (r, z) plane.
        This is primarily used for plotting the source.

        Returns:
            tuple: Tuple (r, z) representing the projection of the source's position in the (r, z) plane.
        """
        x = float(self.position.x)
        y = float(self.position.y)
        z = float(self.position.z)
        
        r = np.sqrt(x**2 + y**2)
        
        return(r, z)
    