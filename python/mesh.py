import os
import pywavefront

class Mesh(object):
    """ Assume only 1 material? (i.e., the texture material?)
    """
    def __init__(self):
        self.path = os.path.normpath("{}{}{}{}".format("Resources", os.sep, name, ".obj"))
        self.faces = []

    def _initialize(self):
        scene = pywavefront.Wavefront(self.path, collect_faces=True)

        # Get material (pywavefront.Wavefront materials contain vertex buffers with texture information)
        # NOTE: We're assuming one material here, so instead of iterating over
        # dict items, we get the first item directly
        mat_name, material = next(iter(scene.materials.items()))

        # Get vertex format string (e.g. 'T2F_N3F_V3F' - tex coords (2 float),
        # normal (3 float), vertex/position (3 float))
        vertex_format = material.vertex_format

        # Get mesh (pywaveFront.Wavefront meshes have faces


