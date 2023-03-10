import sympy


class Layer:
    """
    Class for creating a layer in a tissue model.

    Attributes:
    -----------
    refractive_index: float
        The refractive index of the layer.
    scattering_coefficient: float
        The scattering coefficient of the layer.
    absorption_coefficient: float
        The absorption coefficient of the layer.
    anisotropy_factor: float
        The anisotropy factor of the layer.
    z_range: tuple of two floats
        A tuple of two floats representing the z range of the layer.
    plane: sympy.Plane
        A sympy plane object representing the plane defined by the layer.
    """

    def __init__(self, refractive_index, scattering_coefficient, absorption_coefficient, anisotropy_factor, z_range):
        """
        Constructor for creating a layer.

        Parameters:
        -----------
        refractive_index: float
            The refractive index of the layer.
        scattering_coefficient: float
            The scattering coefficient of the layer.
        absorption_coefficient: float
            The absorption coefficient of the layer.
        anisotropy_factor: float
            The anisotropy factor of the layer.
        z_range: tuple of two floats
            A tuple of two floats representing the z range of the layer.
        """
        self.refractive_index = refractive_index
        self.scattering_coefficient = scattering_coefficient
        self.absorption_coefficient = absorption_coefficient
        self.anisotropy_factor = anisotropy_factor

        self.z_range = z_range
        self.plane = sympy.Plane(sympy.N(sympy.Point3D(
            0, 0, z_range[1])), normal_vector=(0, 0, 1))
