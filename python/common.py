import numpy
import pywavefront
import copy

"""
Notes: 3D vectors are stored (always) as 4-element vectors, for homogeneous coordinates (e.g., translation matrix multiplication)
"""

class Vertex(object):
    def __init__(self, p=None, n=None, t=None):
        # Call the super class constructor (note: the syntax is Python3-specific)
        # In Python2, we would need to call super(Vertex, self).__init__()
        super().__init__()

        # position
        self.p = numpy.array([0,0,0,1]) if p is None else copy.deepcopy(p)

        # normal 
        self.n = numpy.array([0,0,0,1]) if n is None else copy.deepcopy(n)

        # texture coordinates
        self.t = numpy.array([0,0]) if t is None else copy.deepcopy(t)


class Face(object):
    def __init__(self):
        super().__init__()

        self.v0 = Vertex()
        self.v1 = Vertex()
        self.v2 = Vertex()

    def init_by_vertex_copy(self, v0, v1, v2):
        # These assignments assign _REFERENCES_ to existing vertices. TODO: Do we want references, or deep copies? (probably deep copies...)
        self.v0 = copy.deepcopy(v0)
        self.v1 = copy.deepcopy(v1)
        self.v2 = copy.deepcopy(v2)

    def init_by_position(self, p0, p1, p2, face):
        """ Helper function to avoid manually replicating unchanged normals and texture
            coordinates when updating positions
        """
        self.v0 = Vertex(p=p0, n=face.v0.n, t=face.v0.t)
        self.v1 = Vertex(p=p1, n=face.v1.n, t=face.v1.t)
        self.v2 = Vertex(p=p2, n=face.v2.n, t=face.v2.t)

class Color(object):
    """ A class to store normalized colors (color values from 0.0 - 1.0)
    """
    def __init__(self, r=1.0, g=1.0, b=1.0):
        super().__init__()
        self.r = r
        self.g = g
        self.b = b

