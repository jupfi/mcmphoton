import numpy as np 
from timeit import default_timer as timer

class Simulation:
    """
    A Monte Carlo simulation of light transport through tissue.
    
    Attributes:
        type (str): The type of simulation, either "single" or "parallel".
        sim_size (int): The size of the simulation, in millimeters.
        sources (list): The list of sources in the simulation.
        tissue_model (TissueModel): The model of the tissue in the simulation.
    """
    
    def __init__(self, sim_type="single", sim_size = 0.5, multithreaded = False):
        """
        Create a Simulation object.
        
        Args:
            sim_type (str, optional): The type of simulation, either "single" or "package". Default is "single".
            sim_size (int, optional): The radius of the simulation, in centimeters. Default is 10cm.
        """
        self.type = sim_type
        self.sim_size = sim_size
        self.sources = []
        self.tissue_model = None
        self.multithreaded = multithreaded
    
    def add_source(self, source):
        """
        Add a source to the simulation.
        
        Args:
            source (Source): The source to add to the simulation.
        """
        self.sources.append(source)
        
    def add_tissue_model(self, tissue_model):
        """
        Add a tissue model to the simulation. There can only be one tissue model for every simulation.
        
        Args:
            tissue_model (TissueModel): The tissue model to add to the simulation.
        """
        self.tissue_model = tissue_model
    
    # Start the simulation
    def start(self):
        """
        Start the simulation.
        """
        for source in self.sources: 
            source.create_photons()
            source.move_photons_to_boundary()

        
    def simulate(self, timing = False):
        """
        Simulate the light transport through the tissue.
        """
        # Check if there are any photons alive
        start = timer()
        # If multiprocessing is active use the multiprocessing version of the simulation
        if self.multithreaded:
            for source in self.sources:
                source.evaluate_photons()
        else:
            any_alive = np.any(np.array([[photon for photon in source.photons if photon.alive] for source in self.sources]).flatten())
            
            while any_alive:
                for source in self.sources: 
                    source.tick()
                    any_alive = np.any(np.array([[photon for photon in source.photons if photon.alive] for source in self.sources]).flatten())

        end = timer()
        if timing:
            print("The simulation took: {:.2f} seconds".format(end-start))
                
    