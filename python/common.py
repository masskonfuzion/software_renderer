from __future__ import annotations
import numpy

# Notes from numpy documentation (e.g., "help(numpy.ndarray)"):
# |  No ``__init__`` method is needed because the array is fully initialized
# |  after the ``__new__`` method.

# note: the type hinting notation is new to Python as of v3.5
# Today (2019-01-12), it doesn't appear that the types are enforced in any way,
# but the type hinting at least improves readability/understanding, and can be
# used by syntax checkers). See PEP 484
# (https://www.python.org/dev/peps/pep-0484/)

def float_eq( a: float, b: float) -> bool:
    # TODO improve this function -- depending on the inputs, we'll need a
    # "good" way to compute a reasonable EPSILON
    EPSILON = 0.00009

    return ( abs(a - b) <= EPSILON )

def normalize(a: numpy.ndarray) -> numpy.ndarray:
    """ Normalize an array (vector)

        Note that the type hints specify numpy.ndarray -- this is shorthand
        for any of the Point or Vector types in this module (as they all 
        derive from numpy.ndarray
    """
    return a / numpy.linalg.norm(a)

# TODO 2019-02-06 - maybe change the Point2, Point3, Point4, and Matrix4 classes to be lambdas or functions that return numpy.ndarrays, instead of having each class subclass the numpy.ndarray type? Maybe keeps data types cleaner (e.g., right now, if you do m.dot(p), where m is a Matrix4 and p is a Point4, the return type is a Matrix4, when we want Point4; but structurally, the Matrix4 is identical to a Point4 (both are simply numpy.ndarrays, anyway). The disadvantage to this approach is that we lose the .x/y/z/w properties (setters/getters) and the float_eq-modified __eq__ functions -- is there a workaround?


class Point2(numpy.ndarray):
    def __new__(cls, x=0.0, y=0.0):
        """ Initialize the Point2 object

            :param x: X coordinate
            :type x: float
            :param y: Y coordinate
            :type y: float

            according to numpy documentation:
            (e.g., "help(numpy.ndarray)"):

            |  No ``__init__`` method is needed because the array is fully initialized
            |  after the ``__new__`` method.
        """
        return super().__new__(cls, (2,), buffer=numpy.array( (float(x), float(y))  ))

    def __eq__(self, other):
        return float_eq(self[0], other[0]) and \
               float_eq(self[1], other[1])

    ## At one point, we had implemented our own __repr__ -- but we've commented
    ## it out, because our __repr__ was interfering with this class'
    ## compatibility with numpy functions
    ##def __repr__(self):
    ##    return "Point2(x={}, y={})".format(self[0], self[1])

    @property
    def x(self) -> float:
        return self[0]

    @x.setter
    def x(self, value: float):
        self[0] = value

    @property
    def y(self) -> float:
        return self[1]

    @y.setter
    def y(self, value: float):
        self[1] = value

    @staticmethod
    def min(p1, p2):
        """ Return the pariswise minimum of p1 and p2
        """
        return Point2(x=min(p1.x, p2.x), y=min(p1.y, p2.y))

    @staticmethod
    def max(p1, p2):
        """ Return the pariswise maximum of p1 and p2
        """
        return Point2(x=max(p1.x, p2.x), y=max(p1.y, p2.y))

    @staticmethod
    def clamp(p: Point2, lbound: Point2, ubound: Point2) -> Point2:
        """ Clamp the given Point2 to the bounds specified by lbound (lower
            bound) and ubound (upper bound)

            :param p: Point
            :type p: Point2
            :param lbound: Lower boundary
            :type lbound: Point2
            :param ubound: Upper boundary
            :type ubound: Point2
        """
        return Point2( x = max(lbound.x, min(p.x, ubound.x)),
                       y = max(lbound.y, min(p.y, ubound.y)) )


class Point3(numpy.ndarray):
    def __new__(cls, x=0.0, y=0.0, z=0.0):
        """ Initialize the Point3 object
            :param x: X coordinate
            :type x: float
            :param y: Y coordinate
            :type y: float
            :param z: Z coordinate
            :type z: float
    
            # Note that according to numpy documentation:
            |  No ``__init__`` method is needed because the array is fully initialized
            |  after the ``__new__`` method.
        """

        return super().__new__(cls, (3,), buffer=numpy.array( (float(x), float(y), float(z)) ))

    def __eq__(self, other):
        return float_eq(self[0], other[0]) and \
               float_eq(self[1], other[1]) and \
               float_eq(self[2], other[2])

    ## At one point, we had implemented our own __repr__ -- but we've commented
    ## it out, because our __repr__ was interfering with this class'
    ## compatibility with numpy functions
    #def __repr__(self):
    #    return "Point3(x={}, y={}, z={})".format(self[0], self[1], self[2])

    @property
    def x(self) -> float:
        return self[0]

    @x.setter
    def x(self, value: float):
        self[0] = value

    @property
    def y(self) -> float:
        return self[1]

    @y.setter
    def y(self, value: float):
        self[1] = value

    @property
    def z(self) -> float:
        return self[2]

    @z.setter
    def z(self, value: float):
        self[2] = value

    @staticmethod
    def min(p1, p2):
        """ Return the pariswise minimum of p1 and p2
        """
        return Point2(x=min(p1.x, p2.x), y=min(p1.y, p2.y), z=min(p1.z, p2.z))

    @staticmethod
    def max(p1, p2):
        """ Return the pariswise maximum of p1 and p2
        """
        return Point2(x=max(p1.x, p2.x), y=max(p1.y, p2.y), z=max(p1.z, p2.z))


class Point4(numpy.ndarray):
    def __new__(cls, x=0.0, y=0.0, z=0.0, w=0.0):
        """ Initialize the Point4 object

            :param x: X coordinate
            :type x: float
            :param y: Y coordinate
            :type y: float
            :param z: Z coordinate
            :type z: float
            :param w: W coordinate
            :type w: float

            Notes: according to numpy documentation:
            |  No ``__init__`` method is needed because the array is fully initialized
            |  after the ``__new__`` method.
        """
        return super().__new__(cls, (4,), buffer=numpy.array( (float(x), float(y), float(z), float(w)) ))

    def __eq__(self, other):
        return float_eq(self[0], other[0]) and \
               float_eq(self[1], other[1]) and \
               float_eq(self[2], other[2]) and \
               float_eq(self[3], other[3])

    ## At one point, we had implemented our own __repr__ -- but we've commented
    ## it out, because our __repr__ was interfering with this class'
    ## compatibility with numpy functions
    #def __repr__(self):
    #    return "Point4(x={}, y={}, z={}, w={})".format(self[0], self[1], self[2], self[3])

    @property
    def x(self) -> float:
        return self[0]

    @x.setter
    def x(self, value: float):
        self[0] = value

    @property
    def y(self) -> float:
        return self[1]

    @y.setter
    def y(self, value: float):
        self[1] = value

    @property
    def z(self) -> float:
        return self[2]

    @z.setter
    def z(self, value: float):
        self[2] = value

    @property
    def w(self) -> float:
        return self[3]

    @w.setter
    def w(self, value: float):
        self[3] = value

    @staticmethod
    def min(p1, p2):
        """ Return the pariswise minimum of p1 and p2
        """
        return Point2(x=min(p1.x, p2.x), y=min(p1.y, p2.y), z=min(p1.z, p2.z), w=min(p1.w, p2.w))

    @staticmethod
    def max(p1, p2):
        """ Return the pariswise maximum of p1 and p2
        """
        return Point2(x=max(p1.x, p2.x), y=max(p1.y, p2.y), z=max(p1.z, p2.z), w=max(p1.w, p2.w))


class Matrix4(numpy.ndarray):
    """ Matrix class, stored as 4x4 array (column-major)

        i.e. m[0] = [ col0row0 col0row1 col0row2 col0row3 ]
             ...
             m[3] = [ col3row0 col3row1 col3row2 col3row3 ]
    """
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, (4, 4))

    def __init__(self, *args, **kwargs):
        """ Initialize a Matrix4 object
            
            :param args: A tuple of positional parameters
            :type args: float

            This function will process an input (tuple) of 1 or 16 items
            1 item:   Fill the matrix with 16 instances of the 1 value
            16 items: Fill the matrix with the values provided

            Notes: we call super().__init__() by habit, even though we probably
            don't need to. According to numpy documentation (e.g.,
            "help(numpy.ndarray)"):

            |  No ``__init__`` method is needed because the array is fully initialized
            |  after the ``__new__`` method.
        """
        super().__init__()

        if len(args) == 1:
            # Fill the matrix in with args[0]
            for col in range(0, 4):
                for row in range(0, 4):
                    self[col][row] = args[0]
        elif len(args) == 16:
            # IMPORTANT NOTE: numpy stores matrices in row major order. i.e., indices are [row][col]
            # TODO - (re)write __getitem__ functions to allow passing a tuple into the subscript (e.g., self[x,y] instead of self[x][y])
            self[0][0] = args[0]
            self[0][1] = args[1]
            self[0][2] = args[2]
            self[0][3] = args[3]
            
            self[1][0] = args[4]
            self[1][1] = args[5]
            self[1][2] = args[6]
            self[1][3] = args[7]
            
            self[2][0] = args[8]
            self[2][1] = args[9]
            self[2][2] = args[10]
            self[2][3] = args[11]
            
            self[3][0] = args[12]
            self[3][1] = args[13]
            self[3][2] = args[14]
            self[3][3] = args[15]
        elif len(args) not in (1, 16):
            raise Exception("Invalid number of matrix elements provided")
            

    @staticmethod
    def translationMatrix(t: Point3) -> Matrix4:
        """ Return a translation matrix
            
            :param t: The translation vector
            :type t: Point3

            Note: t is a Point3, but is technically a vector (think of t as
            containing the distances to translate, on each of the 3 axes)
        """
        matrix = Matrix4( 1.0, 0.0, 0.0, t.x,
                          0.0, 1.0, 0.0, t.y,
                          0.0, 0.0, 1.0, t.z,
                          0.0, 0.0, 0.0, 1.0 )
        return matrix

    @staticmethod
    def scaleMatrix(s: float) -> Matrix4:
        """ Return a uniform scaling matrix

            :param s: The scale factor
            :type s: float
        """
        matrix = Matrix4( s  , 0.0, 0.0, 0.0,
                          0.0, s  , 0.0, 0.0,
                          0.0, 0.0, s  , 0.0,
                          0.0, 0.0, 0.0, 1.0 )
        return matrix

    @staticmethod
    def rotationMatrix(angle: float, axis: Point3) -> Matrix4:
        """ Return a rotation matrix representing an angle about a given axis
            
            :param angle: The angle in RADIANS
            :type angle: float
            :param axis: The axis about which to rotate
            :type axis: Point3
        """
        u = normalize(axis)
        c = numpy.cos(angle)
        mc = 1.0 - c
        s = numpy.sin(angle)
        xy = u.x * u.y * mc
        xz = u.x * u.z * mc
        yz = u.y * u.z * mc
        xs = u.x * s
        ys = u.y * s
        zs = u.z * s

        matrix = Matrix4( u.x * u.x * mc + c, xy - zs           , xz + ys           , 0.0,
                          xy + zs           , u.y * u.y * mc + c, yz - xs           , 0.0,
                          xz - ys           , yz + xs           , u.z * u.z * mc + c, 0.0,
                          0.0               , 0.0               , 0.0               , 1.0 )

        return matrix

    @staticmethod
    def lookAtMatrix(eye: Point3, target: Point3, up: Point3) -> Matrix4:
        # Note: normalize() returns a "plain" numpy.ndarray; but we want
        # Point3 objects, so we can have .x, .y, and .z properties. So we
        # call the Point3 constructor, passing the result of normalize() as
        # *args, so the constructor can unpack the elements of the ndarray
        n = Point3( *normalize(target - eye) )
        v = normalize(up)
        u = Point3( *normalize(numpy.cross(n, v)) )
        v = Point3( *normalize(numpy.cross(u, n)) )

        dot_u_eye = float( u.dot(eye) )
        dot_v_eye = float( v.dot(eye) )
        dot_n_eye = float( n.dot(eye) )

        matrix = Matrix4(  u.x,  u.y,  u.z, -u.dot(eye),
                           v.x,  v.y,  v.z, -v.dot(eye),
                          -n.x, -n.y, -n.z,  n.dot(eye),
                           0.0,  0.0,  0.0,  1.0 )

        return matrix

    @staticmethod
    def perspectiveMatrix(fov: float, aspect: float, near: float, far: float) -> Matrix4:
        """ Return a perspective matrix
            
            :param fov: Field of view angle in DEGREES
            :type fov: float
            :param aspect: Aspect ratio 
            :type aspect: float
            :param near: Near plane
            :type near: float
            :param far: Far plane
            :type far: float
        """
        matrix = Matrix4(0.0)
        radfov = numpy.pi * fov / 180.0
        f = 1.0 / numpy.tan(radfov / 2.0)

        matrix[0][0] = f / aspect
        matrix[1][1] = f
        matrix[2][2] = (far + near) / (near - far)
        matrix[2][3] = -1.0
        matrix[3][2] = (2.0 * far * near) / (near - far)
        matrix[3][3] = 0.0

        return matrix


class Color(numpy.ndarray):
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, 3)

    def __init__(self, r=0.0, g=0.0, b=0.0):
        """ Initialize the Color object

            This class is basically a fancy container for 3 values

            Notes: we call super().__init__(), even though, according to
            numpy documentation, we might not need to,
            (e.g., "help(numpy.ndarray)"):

            |  No ``__init__`` method is needed because the array is fully initialized
            |  after the ``__new__`` method.
        """
        super().__init__()
        self[0] = float(r)
        self[1] = float(g)
        self[2] = float(b)

    @property
    def r(self) -> float:
        return self[0]

    @r.setter
    def r(self, value: float):
        self[0] = value

    @property
    def g(self) -> float:
        return self[1]

    @g.setter
    def g(self, value: float):
        self[1] = value

    @property
    def b(self) -> float:
        return self[2]

    @b.setter
    def b(self, value: float):
        self[2] = value

    def as256(self):
        """ Return color components as an array with components represented
            in the range 0..255
        """
        return [ int(x * 256) for x in self ]
