import unittest
import numpy as np

from ..classes.point3d import Point3D

class TestPoint3D(unittest.TestCase):
    def setUp(self):
        self.point = Point3D(1, 2, 3)
    
    def test_normalized(self):
        normalized_point = self.point.normalized()
        np.testing.assert_allclose(np.linalg.norm(normalized_point.as_array()), 1)
    
    def test_scaled(self):
        scaled_point = self.point.scaled(2, 3, 4)
        np.testing.assert_allclose(scaled_point.as_array(), np.array([2, 6, 12]))
    
    def test_as_array(self):
        np.testing.assert_allclose(self.point.as_array(), np.array([1, 2, 3]))