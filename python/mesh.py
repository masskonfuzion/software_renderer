
# Import a module that allows annotation (type hinting) of a type within its
# own class definition
from __future__ import annotations
from common import Point2
from common import Point4
from copy import deepcopy
import traceback
import os


class Vertex:
    def __init__(self, p: Point4 = None, n: Point4 = None, t: Point2 = None):
        """ Instantiate a Vertex object

            Note: the data members are all numpy.Point4 objects. The
            PtahRenderer code on which this project is based on uses the
            typealias keyword (in Swift) to easily initialize objects of type,
            Position, Normal, and UV. However, Python does not have a true
            analog for typealias; so instead, I am manually initializing
            numpy.Point4 objects with the appropriate 4D (homogeneous) values
            to represent points and vectors with homogeneous coordinates.
        """
        self.p = deepcopy(p)    # p (position): point (x,y,z,w, where w == 0)
        self.n = deepcopy(n)    # n (normal): vector (x,y,z,w, where w == 1)
        self.t = deepcopy(t)    # t (texture coordinates): point (x,y)

    def __eq__(self, other):
        # Note that this function relies on the == operator of the various
        # Point classes
        return self.p == other.p and \
               self.n == other.n and \
               self.t == other.t

    def __repr__(self):
        return "Vertex(p={},n={},t={})".format(self.p, self.n, self.t)


class Face:
    def __init__(self, v0: Vertex = None, v1: Vertex = None, v2: Vertex = None):
        self.v0 = deepcopy(v0)
        self.v1 = deepcopy(v1)
        self.v2 = deepcopy(v2)

    # Helper function to avoid manually replicating unchanged normals and texture
    # coordinates when updating positions.
    def init(self, p0: Point4, p1: Point4, p2: Point4, face: Face):
        """ Initialize a face

            This function does not have any default parameters
        """
        self.v0 = Vertex(p=p0, n=face.v0.n, t=face.v0.t)
        self.v1 = Vertex(p=p1, n=face.v1.n, t=face.v1.t)
        self.v2 = Vertex(p=p2, n=face.v2.n, t=face.v2.t)

    def __repr__(self):
        return "Face(v0={}, v1={}, v2={})".format(self.v0, self.v1, self.v2)

class Mesh:
    def __init__(self, name):
        self.faces = [] # to be filled in with Face objects

        class FaceIndices:
            def __init__(self, p=None, t=None, n=None):
                """ Initialize a helper object to store the indices (into other
                    arrays) to reference, in order to compose faces in this
                    Mesh 

                    Notes:
                    - In PtahRenderer (written in Swift), this definition is a
                    struct def, not a class. But Python doesn't have a formal
                    struct construct.

                    - In a language with pointers, the entire idea of using
                    FaceIndices might be obviated.
                """
                self.p = p
                self.t = t
                self.n = n

            def __repr__(self):
                return 'FaceIndices(p={}, t={}, n={})'.format(self.p, self.t, self.n)

        path = "{}{}{}{}".format('Resources', os.path.sep, name, '.obj')

        positions = []  # list of Positions (Point4: x,y,z,0.0)
        normals = []    # list of Normals (Point4: x,y,z,1.0)
        uvs = []        # list of UVs (Point2: x,y)
        indices = []    # list of FaceIndices

        try:
            with open(path, 'r') as fd:
                for line in fd:
                    components = line.split()

                    lineprefix = components[0].strip()

                    if lineprefix == "v":
                        # Position (Point4: x,y,z,1)
                        positions.append( Point4(float(components[1].strip()),
                                                 float(components[2].strip()),
                                                 float(components[3].strip()),
                                                 1.0) )
                    elif lineprefix == "vt":
                        # UV coordinates
                        uvs.append( Point2(float(components[1].strip()),
                                           float(components[2].strip())) )

                    elif lineprefix == "vn":
                        # Normal coordinates (Point4: x,y,z,0)
                        normals.append( Point4(float(components[1].strip()),
                                               float(components[2].strip()),
                                               float(components[3].strip())) )

                    elif lineprefix == "f":
                        # Face with positions/uvs/normals. Example:
                        #f  v1/vt1/vn1   v2/vt2/vn2   v3/vt3/vn3 . . .

                        # Rearrange a string, like 
                        #   "f  v1/vt1/vn1   v2/vt2/vn2   v3/vt3/vn3"
                        # into list of lists, like 
                        #   [ [v1, vt1, vn1], [v2, vt2, vn2], [v3, vt3, vn3] ]
                        ids_by_vtx = [s.split("/") for s in line.split(" ")[1:]]
                        #print("Loaded indices by vertex")

                        # Convert the items in the list of lists from string
                        # to int (I have to admit, the Swifty way of doing this
                        # is more concise.. I like it)
                        for i in range(0, len(ids_by_vtx)):
                            # could've said range(len(ids_by_vtx)), (as the
                            # range starts at 0, by default. But the form,
                            # range(0, x) is more explicit
                            ids_by_vtx[i] = list(map(int, ids_by_vtx[i]))
                        #print("Converted indices by vertex (str) into int")

                        ##TODO delete
                        ### - rearrange this "ids_by_vtx" nested list:
                        ###   [ [v1, vt1, vn1], [v2, vt2, vn2], [v3, vt3, vn3] ]
                        ###   into this ids_by_type
                        ###   [ [v1, v2, v3], [vt1, vt2, vt3], [vn1, vn2, vn3] ]
                        ##ids_by_type = list(zip(ids_by_vtx[0], ids_by_vtx[1], ids_by_vtx[2]))
                        ##print("Loaded indices by type")

                        # Construct FaceIndices objects (note: subtract 1 from
                        # each index because .obj files are 1-indexed, but
                        # Python list indices are 0-indexed
                        ###fi = [ FaceIndices(p=ids[0]-1, t=ids[1]-1, n=ids[2]-1) for ids in ids_by_type ] # TODO delete
                        fi = [ FaceIndices(p=ids[0]-1, t=ids[1]-1, n=ids[2]-1) for ids in ids_by_vtx ]
                        indices.append( [fi[0], fi[1], fi[2]] )
                        #print("Loaded face indices")

                        #print("Indices: {}".format(indices))

                # Expand indices, as vertices can have different normals/uvs depending on the face they are linked to.
                self.faces = [ Face( v0=Vertex(p=positions[idx_obj[0].p], n=normals[idx_obj[0].n], t=uvs[idx_obj[0].t]),
                                     v1=Vertex(p=positions[idx_obj[1].p], n=normals[idx_obj[1].n], t=uvs[idx_obj[1].t]),
                                     v2=Vertex(p=positions[idx_obj[2].p], n=normals[idx_obj[2].n], t=uvs[idx_obj[2].t]) ) for idx_obj in indices ]

                #print("Initialized mesh faces")


        except Exception as e:
            print("Couldn't load the mesh: {}".format(e))
            traceback.print_exc()
            return
