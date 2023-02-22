import numpy as np

class Point3D:
    """
    Point3D class represents a point in 3D space with x, y and z coordinates.
    
    Attributes:
        x (float): x coordinate of the point.
        y (float): y coordinate of the point.
        z (float): z coordinate of the point.
    """
    
    def __init__(self, x, y, z):
        """
        Initializes a new Point3D object.
        
        Args:
            x (float): x coordinate of the point.
            y (float): y coordinate of the point.
            z (float): z coordinate of the point.
        """
        self.x = x
        self.y = y
        self.z = z
    
    def normalized(self):
        """
        Normalizes the Point3D object to have length 1.
        
        Returns:
            Point3D: Normalized Point3D object.
        """
        x = self.x / np.sqrt(self.x**2 + self.y**2 + self.z**2)
        y = self.y / np.sqrt(self.x**2 + self.y**2 + self.z**2)
        z = self.z / np.sqrt(self.x**2 + self.y**2 + self.z**2)
        return(Point3D(x,y,z))
    
    def scaled(self, xscale, yscale, zscale):
        """
        Scales the Point3D object by specified amounts in x, y, and z directions.
        
        Args:
            xscale (float): Scale factor for the x coordinate.
            yscale (float): Scale factor for the y coordinate.
            zscale (float): Scale factor for the z coordinate.
        
        Returns:
            Point3D: Scaled Point3D object.
        """
        x = self.x * xscale
        y = self.y * yscale
        z = self.z * zscale
        return(Point3D(x,y,z))
    
    def as_array(self):
        """
        Returns the Point3D object as a numpy array.
        
        Returns:
            numpy.ndarray: 3D numpy array of shape (3,) representing the point.
        """
        return np.array([self.x, self.y, self.z])