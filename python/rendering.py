from common import *
from framebuffer import *
from mesh import *
from math import floor
from texture import *

# Vertex transformation functions

def worldSpaceToClipSpace(faceModelSpace: Face, mvp: Matrix4, lightMatrix: Matrix4) -> Face:
    """ Transform model vertices from world space to clip space
        
        :param faceModelSpace: The face to transform
        :type faceModelSpace: Face
        :param mvp: The MVP (model-view-projection) transformation matrix
        :type mvp: Matrix4
        :param lightMatrix: The lighting transformation matrix
        :type lightMatrix: Matrix4
    """

    # Note the funky syntax here -- it's because I'm using numpy, subclassing
    # numpy's ndarray class for my own purposes
    # numpy's dot() function is how to multiply a matrix by a vector
    # i.e., call it like "mat.dot(v)" ("mat * p" does something different)
    # But... Matrix4.dot(Point4) returns a "Matrix4" (I could write wrapper
    # functions to handle return types better; but because I'm lazy :-D I am
    # constructing a new Point4, passing the result of the dot() function to
    # the constructor as a Python variable-length input (a.k.a. *args)
    v0 = Vertex( p = Point4(*mvp.dot(faceModelSpace.v0.p)),
                 n = Point4(*lightMatrix.dot(faceModelSpace.v0.n)),
                 t = faceModelSpace.v0.t )

    v1 = Vertex( p = Point4(*mvp.dot(faceModelSpace.v1.p)),
                 n = Point4(*lightMatrix.dot(faceModelSpace.v1.n)),
                 t = faceModelSpace.v1.t )

    v2 = Vertex( p = Point4(*mvp.dot(faceModelSpace.v2.p)),
                 n = Point4(*lightMatrix.dot(faceModelSpace.v2.n)),
                 t = faceModelSpace.v2.t )

    return Face( v0=v0, v1=v1, v2=v2 )


def perspectiveDivide(faceClipSpace: Face) -> Face:
    """ Transform model from clip space to normalized space
    """
    w0_inv = 1.0 / faceClipSpace.v0.p.w
    w1_inv = 1.0 / faceClipSpace.v1.p.w
    w2_inv = 1.0 / faceClipSpace.v2.p.w

    p0NormalizedSpace = Point4( faceClipSpace.v0.p.x * w0_inv,
                                faceClipSpace.v0.p.y * w0_inv,
                                faceClipSpace.v0.p.z * w0_inv,
                                faceClipSpace.v0.p.w )

    p1NormalizedSpace = Point4( faceClipSpace.v1.p.x * w1_inv,
                                faceClipSpace.v1.p.y * w1_inv,
                                faceClipSpace.v1.p.z * w1_inv,
                                faceClipSpace.v1.p.w )

    p2NormalizedSpace = Point4( faceClipSpace.v2.p.x * w2_inv,
                                faceClipSpace.v2.p.y * w2_inv,
                                faceClipSpace.v2.p.z * w2_inv,
                                faceClipSpace.v2.p.w )

    # The following code calls init() (a second "constructor"), which is
    # different than the __init__ (constructor) for the Face class. That is
    # because the PtahRenderer's Face class has an overloaded constructor; but
    # Python does not support function overloading
    retface = Face()
    retface.init(p0=p0NormalizedSpace, p1=p1NormalizedSpace, p2=p2NormalizedSpace, face=faceClipSpace)

    return retface

def normalizedSpaceToScreenSpace(faceNormalizedSpace: Face, width: int, height: int) -> Face:
    """ Convert model from normalized space into screen space
    """
    p0ScreenSpace = Point4( floor(0.5 * width * (faceNormalizedSpace.v0.p.x + 1.0)),
                            floor(0.5 * height * (faceNormalizedSpace.v0.p.y + 1.0)),
                            faceNormalizedSpace.v0.p.z,
                            faceNormalizedSpace.v0.p.w )

    p1ScreenSpace = Point4( floor(0.5 * width * (faceNormalizedSpace.v1.p.x + 1.0)),
                            floor(0.5 * height * (faceNormalizedSpace.v1.p.y + 1.0)),
                            faceNormalizedSpace.v1.p.z,
                            faceNormalizedSpace.v1.p.w )

    p2ScreenSpace = Point4( floor(0.5 * width * (faceNormalizedSpace.v2.p.x + 1.0)),
                            floor(0.5 * height * (faceNormalizedSpace.v2.p.y + 1.0)),
                            faceNormalizedSpace.v2.p.z,
                            faceNormalizedSpace.v2.p.w )

    retface = Face()
    retface.init(p0ScreenSpace, p1ScreenSpace, p2ScreenSpace, faceNormalizedSpace)
    return retface


# Drawing functions

def draw(fScreen: Face, framebuffer: Framebuffer, texture: Texture):
    """ Update the frame buffer based on incoming Face data
    """
    # Compute clamped bounding box (mini and maxi are Point2 objs)
    (mini, maxi) = boundingBox(fScreen, framebuffer.width, framebuffer.height)

    for x in range(int(mini.x), int(maxi.x) + 1):
        for y in range(int(mini.y), int(maxi.y) + 1):
            # Compute the barycentric coordinates of the current pixel
            bary = barycentre(Point2(float(x), float(y)), fScreen.v0.p, fScreen.v1.p, fScreen.v2.p)

            # if one of them is negative, we are outside the triangle
            if bary.x < 0.0 or bary.y < 0.0 or bary.z < 0.0:
                continue

            # Interpolate depth at the current pixel
            z = fScreen.v0.p.z * bary.x + fScreen.v1.p.z * bary.y + fScreen.v2.p.z * bary.z

            # If the current triangle pixel is closer than the last one drawn
            if z < framebuffer.getDepth(x, y):
                # Compute perspective-correct interpolation coefficients
                persp = Point3(bary.x/fScreen.v0.p.w, bary.y/fScreen.v1.p.w, bary.z/fScreen.v2.p.w)
                persp = (1.0 / (persp.x + persp.y + persp.z)) * persp

                # Perspective interpolation of texture coordinates and normal
                tex = persp.x * fScreen.v0.t + persp.y * fScreen.v1.t + persp.z * fScreen.v2.t
                nor = persp.x * fScreen.v0.n + persp.y * fScreen.v1.n + persp.z * fScreen.v2.n
                nor = Point3(x=nor.x, y=nor.y, z=nor.z)

                # Compute the color: use a diffuse factor with a strong light intensity
                unit_normal = normalize(nor)
                unit_yz = normalize(Point3(0.0, 1.0, 1.0))

                diffuse = 1.5 * max(0.0, float(unit_normal.transpose().dot(unit_yz)))
                # TODO define texture.. Where does it come from? Is it global, from main?
                #color = diffuse * texture[tex.x][tex.y]
                color = diffuse * texture[tex.x, tex.y]

                # Set color and depth
                framebuffer.setColor(x, y, color)
                framebuffer.setDepth(x, y, z)

        

def boundingBox(vs: Face, width: int, height:int) -> (Point2, Point2):
    """ Return a 2-tuple of Point2 objects, representing a a bounding box around screen space
    """
    # We work only in the (x,y) plane.
    v0ScreenSpace = Point2(vs.v0.p.x, vs.v0.p.y)
    v1ScreenSpace = Point2(vs.v1.p.x, vs.v1.p.y)
    v2ScreenSpace = Point2(vs.v2.p.x, vs.v2.p.y)

    # Find minimal and maximal points.
    # TODO figure out the inner min()/max() calls -- might need to write min/max fns that work on Points
    # PAIRWISE min/max!!!!!
    mini = Point2.min( Point2.min(v0ScreenSpace, v1ScreenSpace), v2ScreenSpace )
    maxi = Point2.max( Point2.max(v0ScreenSpace, v1ScreenSpace), v2ScreenSpace )

    # Framebuffer bounds.
    lim = Point2(float(width) - 1.0, float(height) - 1.0)

    # Clamp the bounding box against the framebuffer.
    finalMin = Point2.clamp(Point2.min(mini, maxi), Point2(0,0), lim)
    finalMax = Point2.clamp(Point2.max(mini, maxi), Point2(0,0), lim)
    return (finalMin, finalMax)

def barycentre(p: Point2, v0: Point4, v1: Point4, v2: Point4) -> Point3:
    # v0 will be the origin.
    # ab and ac: basis vectors
    ab = v1 - v0
    ac = v2 - v0

    # pa: vector with p coordinates in the barycentric frame.
    pa = Point2(v0.x, v0.y) - p

    # This formula for uv1 is pure sorcery
    uv1 = Point3( *numpy.cross(Point3(ac.x, ab.x, pa.x).transpose(), Point3(ac.y, ab.y, pa.y).transpose())[0] )

    # Also note: numpy.cross does not work with 3x1 arrays; it only works
    # (apparently) with 1x3 arrays. My Point3's are 3x1, so I cannot use them
    # in numpy.cross() calculations; I can only use their transpose(). 
    # The result of a numpy.cross, with the data types I've crafted, is a 1x3
    # numpy array
    # Example:
	#>>> vec_a = Point3(1,2,3)
	#>>> vec_b = Point3(4,5,6)
	#>>> numpy.cross(vec_a.transpose(), vec_b.transpose())
	#array([[-3.,  6., -3.]])
    # So the *arg I've passed to the Point3 constructor (above) indexes the
    # [0] item of the numpy.cross() array, which is, itself, an array
    # (I'm not sure why the numpy.cross appears to be actually a 2D array,
    # like [[x, y, z]]. I suppose I could read the documentation, but... Nah)

    # Avoid division imprecision
    if abs(uv1.z) < 1e-2:
        return Point3(-1.0, 1.0, 1.0)

    # Combine
    return (1.0/uv1.z) * Point3( uv1.z - (uv1.x + uv1.y), uv1.y, uv1.x )
