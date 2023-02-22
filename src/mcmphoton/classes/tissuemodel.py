import math

from .layer import Layer

class TissueModel:
    """
    This class represents a tissue model, which is a collection of layered media with specified optical properties.
    
    Attributes:
        layers (List[Layer]): A list of Layer objects representing the layered media in the tissue model.
    """
    def __init__(self):
        """
        Initializes the TissueModel object. The first layer is air and extends to '-infinity'
        """
        first_layer = Layer(1, 1, 1, 0, (-math.inf, 0))
        
        self.layers = [first_layer]
    
    #All other layers stacked between the first_layer and the last_layer
    def add_layer(self, layer):
        """
        Adds a layer to the tissue model and appends a new last layer of air.
        
        Args:
            layer (Layer): The Layer object to be added to the tissue model.
        """
        if len(self.layers) > 1:
            self.layers.pop()
        
        self.layers.append(layer)
        boundary = layer.z_range[1]
        last_layer =Layer(1, 1, 1, 0, (boundary, boundary + 5))
        self.layers.append(last_layer)
        