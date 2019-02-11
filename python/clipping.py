
from mesh import Face, Vertex

def clip(faceClipSpace: Face) -> [Face]:

    # Literally, I didn't know Python supported type hinting like this
    # I found out by doing this project :-)
    triangles : [Face] = []

    # All vertices are behind the camera
    if faceClipSpace.v0.p.w <= 0.0 and faceClipSpace.v1.p.w <= 0.0 and faceClipSpace.v2.p.w <= 0.0:
        return []

    # All vertices are in front of the camera and inside the frustum
    if faceClipSpace.v0.p.w > 0.0 and faceClipSpace.v1.p.w > 0.0 and faceClipSpace.v2.p.w > 0.0 and \
       abs(faceClipSpace.v0.p.z) < faceClipSpace.v0.p.w and \
       abs(faceClipSpace.v1.p.z) < faceClipSpace.v1.p.w and \
       abs(faceClipSpace.v2.p.z) < faceClipSpace.v2.p.w:
        triangles = [ faceClipSpace ]   # a list with a single face in it

    else:
        # Clip each edge, accumulating vertices that we add or keep in an array
        vertices : [Vertex] = []

        # Note that Python treats a list passed into a function as a pass-by-reference
        # (the PtahRenderer code (in Swift, pass-by-reference must be explicitly declared)
        clipEdge(v0=faceClipSpace.v0, v1=faceClipSpace.v1, vertices=vertices)
        clipEdge(v0=faceClipSpace.v1, v1=faceClipSpace.v2, vertices=vertices)
        clipEdge(v0=faceClipSpace.v2, v1=faceClipSpace.v0, vertices=vertices)

        # If not enough vertices to create a triangular face
        if len(vertices) < 3:
            return []

        # We potentially have a duplicate at the end, that we can remove
        if vertices[len(vertices) - 1] == vertices[0]:
            vertices.pop()

        # Generate a fan of triangles, all sharing the first vertex
        for i in range(1, len(vertices)-1):
            triangles.append(Face(v0=vertices[0], v1=vertices[i], v2=vertices[i+1]))

    return triangles


def clipEdge(v0: Vertex, v1: Vertex, vertices: [Vertex]):
    v0New = v0
    v1New = v1

    # Testing with respect to near plane
    v0Inside = v0.p.w > 0.0 and v0.p.z > -v0.p.w
    v1Inside = v1.p.w > 0.0 and v1.p.z > -v1.p.w

    if v0Inside and v1Inside:
        ## Great, nothing to do :-D
        pass
    elif v0Inside or v1Inside:
        # Compute interpolation coefficients
        d0 = v0.p.z + v0.p.w
        d1 = v1.p.z + v1.p.w
        factor = 1.0 / (d1 - d0)
        # New vertex with interpolated coefficients
        newVertex = Vertex( p=factor * (d1 * v0.p - d0 * v1.p),
                            n=factor * (d1 * v0.n - d0 * v1.n),
                            t=factor * (d1 * v0.t - d0 * v1.t) )

        if v0Inside:
            v1New = newVertex
        else:
            v0New = newVertex
    else:
        # Both are outside, on the same side. Remove the edge
        return

    # Add the first vertex if not already added
    num_vertices = len(vertices)
    if num_vertices == 0 or not (vertices[num_vertices - 1] == v0New):
        # "not (v1 == v2)", instead of "v1 != v2", because we never
        # implemented a "!=" function.. I know... lazy :-)
        vertices.append(v0New)

    # Add the second vertex
    vertices.append(v1New)


def cullFace(faceNormalizedSpace: Face) -> bool:
    d = (faceNormalizedSpace.v1.p.x - faceNormalizedSpace.v0.p.x) * \
        (faceNormalizedSpace.v2.p.y - faceNormalizedSpace.v0.p.y) - \
        (faceNormalizedSpace.v1.p.y - faceNormalizedSpace.v0.p.y) * \
        (faceNormalizedSpace.v2.p.x - faceNormalizedSpace.v0.p.x)
    return d < 0.0
