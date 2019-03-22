#pragma once
#include <iostream>
#include <array>
#include <cmath>
#include <string>
#include <vector>

// Without using C++17 shenanigans, I'm not sure of a better way to assign
// the os path separator. (Realistically, anyone who compiles this code
// probably has access to a C++17 compiler, but also.. can't be too sure)
const std::string os_path_sep =
#ifdef _WIN32
    "\\";
#else
    "/";
#endif


// Split the string, s, on char c; store the resulting list of substrings in v_out
// NOTE: this function is destructive - v_out will always be cleared out
// NOTE ALSO: this function has not been tested on strings with consecutive split_char characters
void string_split(const std::string& s, const char& split_char, std::vector<std::string>& v_out) {
    // Or we could make this function return an int "success" code
    v_out.clear();

    std::size_t startidx = 0;
    //std::cout << "split_string(\"" << s << "\", '" << split_char << "') set startidx: " << startidx << std::endl;
    std::size_t split_char_idx = std::string::npos;

    while (startidx != std::string::npos) {
        split_char_idx = s.find(split_char, startidx);
        //std::cout << "split_string(\"" << s << "\", '" << split_char << "') found split_char_idx: " << split_char_idx << std::endl;
        if (split_char_idx == std::string::npos) {
            split_char_idx = s.length();
            //std::cout << "split_string(\"" << s << "\", '" << split_char << "') adjusted split_char_idx: " << split_char_idx << std::endl;
        }
        
        // Boundary case check -- if split_char_idx == startidx, they both must be 0, which
        // means the first char of the string is the split_char. In that case, the first
        // substr in the vector would be "". But, don't want/need an empty string as the first item
        // in the output vector, so we prevent it from being inserted into the vector.
        if (split_char_idx != startidx) {   
            //std::cout << "split_string(\"" << s << "\", '" << split_char << "') found substr: " << s.substr(startidx, split_char_idx - startidx) << std::endl;
            v_out.push_back( s.substr(startidx, split_char_idx - startidx)  );
        }

        startidx = split_char_idx + 1;
        //std::cout << "split_string(\"" << s << "\", '" << split_char << "') set startidx: " << startidx << std::endl;
        if (startidx > s.length() - 1) {
            startidx = std::string::npos;
            //std::cout << "split_string(\"" << s << "\", '" << split_char << "') adjusted startidx: " << startidx << std::endl;
        }
    }
}


// Type aliases (easier-to-read than typedef)

bool float_eq( float a, float b) {
    // TODO improve this function -- depending on the inputs, we'll need a
    // "good" way to compute a reasonable EPSILON
    float EPSILON = 0.00009;

    return std::abs(a - b) <= EPSILON;
}

// free subtraction operator for Point2, Point3, and Point4
// NOTE/TODO - possibly change the free subtraction operator into a class member?
template <class T>
T operator- (const T& a, const T& b) {
    T ret;
    for (int i = 0; i < ret.v.size(); i++) {
        ret.v[i] = a.v[i] - b.v[i];
    }

    return ret;
}

// free addition operator for Point2, Point3, and Point4
// NOTE/TODO - possibly change the free subtraction operator into a class member?
template <class T>
T operator+ (const T& a, const T& b) {
    T ret;
    for (int i = 0; i < ret.v.size(); i++) {
        ret.v[i] = a.v[i] + b.v[i];
    }

    return ret;
}


// Templated dot product
template <class T>
float dot(const T& a, const T& b) {
    float ret = 0.0;
    for (int i = 0; i < a.v.size(); i++) {
       ret += a[i] * b[i]; 
    }

    return ret;
}


// A templated normalize function
// Normalize a Point2, Point3 or Point4
template <class T>
T normalize(const T& vec) {
    // A templated function to normalize vectors
    // Hopefully the compiler can optimize this into crazy fast ish (SIMD)
    float invlen = float(1.0) / veclen(vec);
    T retvec;
    for (int i = 0; i < retvec.v.size(); i++) {
        retvec[i] = vec[i] * invlen;
    }
    return retvec;
}

template <class T>
float veclen(const T& vec) {
    // A tepmlated function to calculate the length of a vector
    // Uses a loop. Hopefully, the compiler can optimize this heavily (SIMD)
    
    float sqrlen = 0.0;

    for (int i = 0; i < vec.v.size(); i++) {
        sqrlen += vec[i] * vec[i];
    }

    return std::sqrt(sqrlen);
}


struct Color {
    Color(float r = 0.0, float g = 0.0, float b = 0.0) {
        v[0] = r;
        v[1] = g;
        v[2] = b;
    }

    Color(const Color& src) {
        // I should simply be able to copy the other color's items
        for (int i = 0; i < src.v.size(); i++) {
            v[i] = src.v[i];
        }

        // Initialize references to array items;
        r = v[0];
        g = v[1];
        b = v[2];
    }

    // Member functions
    //// TODO see if you can initialize references, such that, for example, myColor.g is a reference to myColor.v[1]
    //float getr() const { return v[0]; }
    //void setr(float r) { v[0] = r; }

    //float getg() const { return v[1]; }
    //void setg(float g) { v[1] = g; }

    //float getb() const { return v[2]; }
    //void setb(float b) { v[2] = b; }

    // Copy assignment operator
    Color& operator= (const Color& color) {
        v[0] = color.v[0];
        v[1] = color.v[1];
        v[2] = color.v[2];

        return *this;
    }

    // TODO delete -- this operator is overloaded by a free operator (also -- document that operator)
    //// Return a Color, the result of multiplying this Color by a scalar
    //Color operator* (float scalar) {
    //    return Color(scalar*v[0], scalar*v[1], scalar*v[2]);
    //}

    std::array<float, 3> as256() {
        return std::array<float, 3> { v[0]*256, v[1]*256, v[2]*256 };
    }


    // the array that contains the color's r, g, and b components
    std::array<float, 3> v;

    float& r = v[0];
    float& g = v[1];
    float& b = v[2];

};


// TODO delete below
//// Color * float scalar multiplication
//// I had to write a separate multiplication function, because I named the
//// underlying array c, instead of what I used for Points (v)
//Color operator* (float scalar, const Color& vec) {
//    Color ret;
//    for (int i = 0; i < vec.v.size(); i++) {
//        ret.v[i] = vec.v[i] * scalar;
//    }
//    return ret;
//}
//
//template <class Color>
//Color operator* (const Color& vec, float scalar) {
//    Color ret;
//    for (int i = 0; i < vec.v.size(); i++) {
//        ret.v[i] = vec.v[i] * scalar;
//    }
//    return ret;
//}


// A 2D point/vector data type
// NOTE: this struct uses scalar data types.
// TODO come back later and replace with SIMD-ready vectorized types
struct Point2 {
    Point2(float x = 0.0, float y = 0.0) {
        v[0] = x;
        v[1] = y;
    }

    Point2(const Point2& src) {
        for (int i = 0; i < src.v.size(); i++) {
            v[i] = src.v[i];
        }

        // Initialize references to array items;
        x = v[0];
        y = v[1];
    }

    // Member functions
    float getx() const { return v[0]; }
    void setx(float x) { v[0] = x; }

    float gety() const { return v[1]; }
    void sety(float y) { v[1] = y; }

    // Accessor function - bounds checking will (hopefully) be handled by the std::array class
    const float& operator[] (int i) const { return v[i]; }

    // Mutator function - bounds checking will (hopefully) be handled by the std::array class
    float& operator[] (int i) { return v[i]; }

    // Copy assignment operator
    Point2& operator= (const Point2& point) {
        v[0] = point.v[0];
        v[1] = point.v[1];
        return *this;
    }

    static Point2 min(const Point2& p1, const Point2& p2) {
        return Point2( std::min(p1.v[0], p2.v[0]), std::min(p1.v[1], p2.v[1]) );
    }

    static Point2 max(const Point2& p1, const Point2& p2) {
        return Point2( std::max(p1.v[0], p2.v[0]), std::max(p1.v[1], p2.v[1]) );
    }

    static Point2 clamp(Point2 p, Point2 lbound, Point2 ubound) {
        //Clamp the given Point2 to the bounds specified by lbound (lower
        //bound) and ubound (upper bound)

        //:param p: Point
        //:param lbound: Lower boundary
        //:param ubound: Upper boundary
        
        return Point2( std::max(lbound.getx(), std::min(p.getx(), ubound.getx())),
                       std::max(lbound.gety(), std::min(p.gety(), ubound.gety())) );
        
    }

    // Data members
    std::array< float, 2> v;

    float& x = v[0];
    float& y = v[1];
};


// A 3D point/vector data type
// NOTE: this struct uses scalar data types.
// TODO come back later and replace with SIMD-ready vectorized types
struct Point3 {
    Point3(float x = 0.0, float y = 0.0, float z = 0.0) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    Point3(const Point3& src) {
        for (int i = 0; i < src.v.size(); i++) {
            v[i] = src.v[i];
        }

        // Initialize references to array items;
        x = v[0];
        y = v[1];
        z = v[2];
        
    }

    // Member functions
    float getx() const { return v[0]; }
    void setx(float x) { v[0] = x; }

    float gety() const { return v[1]; }
    void sety(float y) { v[1] = y; }

    float getz() const { return v[2]; }
    void setz(float z) { v[2] = z; }

    // Accessor function - bounds checking will (hopefully) be handled by the std::array class
    const float& operator[] (int i) const { return v[i]; }

    // Mutator function - bounds checking will (hopefully) be handled by the std::array class
    float& operator[] (int i) { return v[i]; }

    // Copy assignment operator
    Point3& operator= (const Point3& point) {
        v[0] = point.v[0];
        v[1] = point.v[1];
        v[2] = point.v[2];
        return *this;
    }

    // Data members
    std::array< float, 3> v;

    float& x = v[0];
    float& y = v[1];
    float& z = v[2];
};

// Cross product
Point3 cross(const Point3& a, const Point3& b) {
    // Return a Point3 
    return Point3( a.v[1] * b.v[2] - a.v[2] * b.v[1],
    			   a.v[2] * b.v[0] - a.v[0] * b.v[2],
    			   a.v[0] * b.v[1] - a.v[1] * b.v[0] );
}


// A 4D point/vector data type (3D + homogeneous coordinate)
// NOTE: this struct uses scalar data types.
// TODO come back later and replace with SIMD-ready vectorized types
struct Point4 {
    Point4(float x = 0.0, float y = 0.0, float z = 0.0, float w = 0.0) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
        v[3] = w;
    }

    Point4(const Point4& src) {
        for (int i = 0; i < src.v.size(); i++) {
            v[i] = src.v[i];
        }

        // Initialize references to array items;
        x = v[0];
        y = v[1];
        z = v[2];
        w = v[3];

    }

    // Member functions
    float getx() const { return v[0]; }
    void setx(float x) { v[0] = x; }

    float gety() const { return v[1]; }
    void sety(float y) { v[1] = y; }

    float getz() const { return v[2]; }
    void setz(float z) { v[2] = z; }

    float getw() const { return v[3]; }
    void setw(float w) { v[3] = w; }

    // Accessor function - bounds checking will (hopefully) be handled by the std::array class
    const float& operator[] (int i) const { return v[i]; }

    // Mutator function - bounds checking will (hopefully) be handled by the std::array class
    float& operator[] (int i) { return v[i]; }

    // Copy assignment operator
    Point4& operator= (const Point4& point) {
        v[0] = point.v[0];
        v[1] = point.v[1];
        v[2] = point.v[2];
        v[3] = point.v[3];
        return *this;
    }

    // Data members
    std::array< float, 4> v;

    float& x = v[0];
    float& y = v[1];
    float& z = v[2];
    float& w = v[3];
};




// A 4x4 matrix data type
// NOTE: this struct uses scalar data types.
// TODO come back later and replace with SIMD-ready vectorized types
struct Matrix4 {
    // NOTE: Matrix4 class is column major
    // i.e., this->m is a 4x4 array -- it is an array of 4 "things"; each
	// "thing" is a 4D vector, which represents a column of the matrix
    // As a result, vectors will be post-multiplied (v' = Mv)
    Matrix4() {
        // Default constructor initializes the identity matrix
        // column 0
        m[0][0] = 1.0;
        m[0][1] = 0.0;
        m[0][2] = 0.0;
        m[0][3] = 0.0;

        // column 1
        m[1][0] = 0.0;
        m[1][1] = 1.0;
        m[1][2] = 0.0;
        m[1][3] = 0.0;

        // column 2
        m[2][0] = 0.0;
        m[2][1] = 0.0;
        m[2][2] = 1.0;
        m[2][3] = 0.0;

        // column 3
        m[3][0] = 0.0;
        m[3][1] = 0.0;
        m[3][2] = 0.0;
        m[3][3] = 1.0;
    }

    Matrix4 (const Matrix4& src) {
        // Explicitly call the copy constructors for each Point4 in the m array
        m[0] = Point4(src.m[0]);
        m[1] = Point4(src.m[1]);
        m[2] = Point4(src.m[2]);
        m[3] = Point4(src.m[3]);
    }

    Matrix4(float f) {
        // Initialize all values to f
        // column 0
        m[0][0] = f;
        m[0][1] = f;
        m[0][2] = f;
        m[0][3] = f;

        // column 1
        m[1][0] = f;
        m[1][1] = f;
        m[1][2] = f;
        m[1][3] = f;

        // column 2
        m[2][0] = f;
        m[2][1] = f;
        m[2][2] = f;
        m[2][3] = f;

        // column 3
        m[3][0] = f;
        m[3][1] = f;
        m[3][2] = f;
        m[3][3] = f;
    }

    Matrix4(float m00, float m10, float m20, float m30,
            float m01, float m11, float m21, float m31,
            float m02, float m12, float m22, float m32,
            float m03, float m13, float m23, float m33) {

        // Specify the values the way they'd look if you wrote the matrix on
        // paper (the way they teach it in school)
        
        // Indices are m[column][index_within_column]
        // We write the matrix the way it looks on paper, but under the hood,
        // we store each column as a contiguous array
        m[0][0] = m00;  m[1][0] = m10;  m[2][0] = m20;  m[3][0] = m30;
        m[0][1] = m01;  m[1][1] = m11;  m[2][1] = m21;  m[3][1] = m31;
        m[0][2] = m02;  m[1][2] = m12;  m[2][2] = m22;  m[3][2] = m32;
        m[0][3] = m03;  m[1][3] = m13;  m[2][3] = m23;  m[3][3] = m33;
    }


    // Class member functions
    // Accessor function - bounds checking will (hopefully) be handled by the std::array class
    const Point4& operator[] (int i) const { return m[i]; }

    // Mutator function - bounds checking will (hopefully) be handled by the std::array class
    Point4& operator[] (int i) { return m[i]; }

    static Matrix4 translationMatrix(float x, float y, float z) {
        return Matrix4( 1.0, 0.0, 0.0, x,
                        0.0, 1.0, 0.0, y,
                        0.0, 0.0, 1.0, z,
                        0.0, 0.0, 0.0, 1.0 );
    }

    static Matrix4 translationMatrix(const Point3& t) {
        // Overloaded
        return translationMatrix(t[0], t[1], t[2]);
    }

    static Matrix4 scaleMatrix(float s) {
        // Uniform scale
        return Matrix4( s  , 0.0, 0.0, 0.0,
                        0.0, s  , 0.0, 0.0,
                        0.0, 0.0, s  , 0.0,
                        0.0, 0.0, 0.0, 1.0 );
    }

    static Matrix4 rotationMatrix(float angle, Point3 axis) {
        //Return a rotation matrix representing an angle about a given axis
        //    
        //    :param angle: The angle in RADIANS
        //    :param axis: The axis about which to rotate

        Point3 u = normalize(axis);
        float c = std::cos(angle);
        float mc = static_cast<float>(1.0) - c;
        float s = std::sin(angle);
        float xy = u.getx() * u.gety() * mc;
        float xz = u.getx() * u.getz() * mc;
        float yz = u.gety() * u.getz() * mc;
        float xs = u.getx() * s;
        float ys = u.gety() * s;
        float zs = u.getz() * s;

        return Matrix4( u.getx() * u.getx() * mc + c, xy - zs                     , xz + ys                     , 0.0,
                        xy + zs                     , u.gety() * u.gety() * mc + c, yz - xs                     , 0.0,
                        xz - ys                     , yz + xs                     , u.getz() * u.getz() * mc + c, 0.0,
                        0.0                         , 0.0                         , 0.0                         , 1.0 );
    }

    static Matrix4 lookAtMatrix(Point3 eye, Point3 target, Point3 up) {
        Point3 n = normalize(target - eye); // TODO...we might need to write this. Our Point and Vector classes might need to contain the SIMD-friendly data types we want to work with; but this might not work if those classes ARE the raw SIMD-friendly classes
        Point3 v = normalize(up);
        Point3 u = normalize(cross(n, v));    // TODO write cross product (free operator -- no need for a friend function because we're working with a struct -- public by default)
        v = normalize(cross(u, n));

        return Matrix4 (  u.getx(),  u.gety(),  u.getz(), -dot(u, eye),
                          v.getx(),  v.gety(),  v.getz(), -dot(v, eye),
                         -n.getx(), -n.gety(), -n.getz(),  dot(n, eye),
                          0.0     ,  0.0     ,  0.0     ,  1.0 );
    }

    static Matrix4 perspectiveMatrix(float fov, float aspect, float near, float far) {
        //Return a perspective matrix
        //
        //:param fov: Field of view angle in DEGREES
        //:param aspect: Aspect ratio 
        //:param near: Near plane
        //:param far: Far plane
        Matrix4 matrix(0.0);

        float radfov = M_PI * fov / 180.0;
        float f = 1.0 / std::tan(radfov / 2.0);
        
        matrix[0][0] = f / aspect;
        matrix[1][1] = f;
        matrix[2][2] = (far + near) / (near - far);
        matrix[2][3] = (2.0 * far * near) / (near - far);
        matrix[3][2] = -1.0;
        matrix[3][3] = 0.0;

        return matrix;
    }

    // Return the inverse of the given matrix (if the given matrix is
    // invertible); else, return a 0 matrix
    static Matrix4 inverse(const Matrix4& m) {
        Matrix4 ret;
        float det;  // determinant

        // m[0][0]..[3] map to mat_col_maj[0]..[3]
        // m[1][0]..[3] map to mat_col_maj[4]..[7]
        // m[2][0]..[3] map to mat_col_maj[8]..[11]
        // m[3][0]..[3] map to mat_col_maj[12]..[15]

        ret[0][0] = m[1][1]  * m[2][2] * m[3][3] - 
                    m[1][1]  * m[2][3] * m[3][2] - 
                    m[2][1]  * m[1][2]  * m[3][3] + 
                    m[2][1]  * m[1][3]  * m[3][2] +
                    m[3][1] * m[1][2]  * m[2][3] - 
                    m[3][1] * m[1][3]  * m[2][2];

        ret[1][0] = -m[1][0]  * m[2][2] * m[3][3] + 
                     m[1][0]  * m[2][3] * m[3][2] + 
                     m[2][0]  * m[1][2]  * m[3][3] - 
                     m[2][0]  * m[1][3]  * m[3][2] - 
                     m[3][0] * m[1][2]  * m[2][3] + 
                     m[3][0] * m[1][3]  * m[2][2];

        ret[2][0] = m[1][0]  * m[2][1] * m[3][3] - 
                    m[1][0]  * m[2][3] * m[3][1] - 
                    m[2][0]  * m[1][1] * m[3][3] + 
                    m[2][0]  * m[1][3] * m[3][1] + 
                    m[3][0] * m[1][1] * m[2][3] - 
                    m[3][0] * m[1][3] * m[2][1];

        ret[3][0] = -m[1][0]  * m[2][1] * m[3][2] + 
                     m[1][0]  * m[2][2] * m[3][1] +
                     m[2][0]  * m[1][1] * m[3][2] - 
                     m[2][0]  * m[1][2] * m[3][1] - 
                     m[3][0] * m[1][1] * m[2][2] + 
                     m[3][0] * m[1][2] * m[2][1];

        ret[0][1] = -m[0][1]  * m[2][2] * m[3][3] + 
                     m[0][1]  * m[2][3] * m[3][2] + 
                     m[2][1]  * m[0][2] * m[3][3] - 
                     m[2][1]  * m[0][3] * m[3][2] - 
                     m[3][1] * m[0][2] * m[2][3] + 
                     m[3][1] * m[0][3] * m[2][2];

        ret[1][1] = m[0][0]  * m[2][2] * m[3][3] - 
                    m[0][0]  * m[2][3] * m[3][2] - 
                    m[2][0]  * m[0][2] * m[3][3] + 
                    m[2][0]  * m[0][3] * m[3][2] + 
                    m[3][0] * m[0][2] * m[2][3] - 
                    m[3][0] * m[0][3] * m[2][2];

        ret[2][1] = -m[0][0]  * m[2][1] * m[3][3] + 
                     m[0][0]  * m[2][3] * m[3][1] + 
                     m[2][0]  * m[0][1] * m[3][3] - 
                     m[2][0]  * m[0][3] * m[3][1] - 
                     m[3][0] * m[0][1] * m[2][3] + 
                     m[3][0] * m[0][3] * m[2][1];

        ret[3][1] = m[0][0]  * m[2][1] * m[3][2] - 
                    m[0][0]  * m[2][2] * m[3][1] - 
                    m[2][0]  * m[0][1] * m[3][2] + 
                    m[2][0]  * m[0][2] * m[3][1] + 
                    m[3][0] * m[0][1] * m[2][2] - 
                    m[3][0] * m[0][2] * m[2][1];

        ret[0][2] = m[0][1]  * m[1][2] * m[3][3] - 
                    m[0][1]  * m[1][3] * m[3][2] - 
                    m[1][1]  * m[0][2] * m[3][3] + 
                    m[1][1]  * m[0][3] * m[3][2] + 
                    m[3][1] * m[0][2] * m[1][3] - 
                    m[3][1] * m[0][3] * m[1][2];

        ret[1][2] = -m[0][0]  * m[1][2] * m[3][3] + 
                     m[0][0]  * m[1][3] * m[3][2] + 
                     m[1][0]  * m[0][2] * m[3][3] - 
                     m[1][0]  * m[0][3] * m[3][2] - 
                     m[3][0] * m[0][2] * m[1][3] + 
                     m[3][0] * m[0][3] * m[1][2];

        ret[2][2] = m[0][0]  * m[1][1] * m[3][3] - 
                    m[0][0]  * m[1][3] * m[3][1] - 
                    m[1][0]  * m[0][1] * m[3][3] + 
                    m[1][0]  * m[0][3] * m[3][1] + 
                    m[3][0] * m[0][1] * m[1][3] - 
                    m[3][0] * m[0][3] * m[1][1];

        ret[3][2] = -m[0][0]  * m[1][1] * m[3][2] + 
                     m[0][0]  * m[1][2] * m[3][1] + 
                     m[1][0]  * m[0][1] * m[3][2] - 
                     m[1][0]  * m[0][2] * m[3][1] - 
                     m[3][0] * m[0][1] * m[1][2] + 
                     m[3][0] * m[0][2] * m[1][1];

        ret[0][3] = -m[0][1] * m[1][2] * m[2][3] + 
                     m[0][1] * m[1][3] * m[2][2] + 
                     m[1][1] * m[0][2] * m[2][3] - 
                     m[1][1] * m[0][3] * m[2][2] - 
                     m[2][1] * m[0][2] * m[1][3] + 
                     m[2][1] * m[0][3] * m[1][2];

        ret[1][3] = m[0][0] * m[1][2] * m[2][3] - 
                    m[0][0] * m[1][3] * m[2][2] - 
                    m[1][0] * m[0][2] * m[2][3] + 
                    m[1][0] * m[0][3] * m[2][2] + 
                    m[2][0] * m[0][2] * m[1][3] - 
                    m[2][0] * m[0][3] * m[1][2];

        ret[2][3] = -m[0][0] * m[1][1] * m[2][3] + 
                     m[0][0] * m[1][3] * m[2][1] + 
                     m[1][0] * m[0][1] * m[2][3] - 
                     m[1][0] * m[0][3] * m[2][1] - 
                     m[2][0] * m[0][1] * m[1][3] + 
                     m[2][0] * m[0][3] * m[1][1];

        ret[3][3] = m[0][0] * m[1][1] * m[2][2] - 
                    m[0][0] * m[1][2] * m[2][1] - 
                    m[1][0] * m[0][1] * m[2][2] + 
                    m[1][0] * m[0][2] * m[2][1] + 
                    m[2][0] * m[0][1] * m[1][2] - 
                    m[2][0] * m[0][2] * m[1][1];

        det = m[0][0] * ret[0][0] + m[0][1] * ret[1][0] + m[0][2] * ret[2][0] + m[0][3] * ret[3][0];

        // TODO replace with float/epsilon comparison
        if (det == 0)
            // If det == 0, then the inversion failed. Return a 0 matrix
            return Matrix4(0.0);

        det = 1.0 / det;

        for (int c = 0; c < 4; c++) {
            for (int r = 0; r < 4; r++) {
                ret[c][r] *= det;
            }
        }

        return ret;
    }

    static Matrix4 transpose(const Matrix4& m) {
        // Store our columns as rows in the new Matrix, ret
        return Matrix4 (m[0][0], m[0][1], m[0][2], m[0][3],
                        m[1][0], m[1][1], m[1][2], m[1][3],
                        m[2][0], m[2][1], m[2][2], m[2][3],
                        m[3][0], m[3][1], m[3][2], m[3][3]);
    }

    // Member variables
    //std::array< std::array<float, 4>, 4 > m
    std::array<Point4, 4> m;

};


// Matrix4 * Matrix4 multiplication
Matrix4 operator* (const Matrix4& a, const Matrix4& b) {
    Matrix4 ret;

    // row 1
    ret[0][0] = a[0][0]*b[0][0] + a[1][0]*b[0][1] + a[2][0]*b[0][2] + a[3][0]*b[0][3];
    ret[1][0] = a[0][0]*b[1][0] + a[1][0]*b[1][1] + a[2][0]*b[1][2] + a[3][0]*b[1][3];
    ret[2][0] = a[0][0]*b[2][0] + a[1][0]*b[2][1] + a[2][0]*b[2][2] + a[3][0]*b[2][3];
    ret[3][0] = a[0][0]*b[3][0] + a[1][0]*b[3][1] + a[2][0]*b[3][2] + a[3][0]*b[3][3];
    
    // row 2
    ret[0][1] = a[0][1]*b[0][0] + a[1][1]*b[0][1] + a[2][1]*b[0][2] + a[3][1]*b[0][3];
    ret[1][1] = a[0][1]*b[1][0] + a[1][1]*b[1][1] + a[2][1]*b[1][2] + a[3][1]*b[1][3];
    ret[2][1] = a[0][1]*b[2][0] + a[1][1]*b[2][1] + a[2][1]*b[2][2] + a[3][1]*b[2][3];
    ret[3][1] = a[0][1]*b[3][0] + a[1][1]*b[3][1] + a[2][1]*b[3][2] + a[3][1]*b[3][3];
    
    // row 3
    ret[0][2] = a[0][2]*b[0][0] + a[1][2]*b[0][1] + a[2][2]*b[0][2] + a[3][2]*b[0][3];
    ret[1][2] = a[0][2]*b[1][0] + a[1][2]*b[1][1] + a[2][2]*b[1][2] + a[3][2]*b[1][3];
    ret[2][2] = a[0][2]*b[2][0] + a[1][2]*b[2][1] + a[2][2]*b[2][2] + a[3][2]*b[2][3];
    ret[3][2] = a[0][2]*b[3][0] + a[1][2]*b[3][1] + a[2][2]*b[3][2] + a[3][2]*b[3][3];
    
    // row 4
    ret[0][3] = a[0][3]*b[0][0] + a[1][3]*b[0][1] + a[2][3]*b[0][2] + a[3][3]*b[0][3];
    ret[1][3] = a[0][3]*b[1][0] + a[1][3]*b[1][1] + a[2][3]*b[1][2] + a[3][3]*b[1][3];
    ret[2][3] = a[0][3]*b[2][0] + a[1][3]*b[2][1] + a[2][3]*b[2][2] + a[3][3]*b[2][3];
    ret[3][3] = a[0][3]*b[3][0] + a[1][3]*b[3][1] + a[2][3]*b[3][2] + a[3][3]*b[3][3];

    return ret;
};


// Matrix4 * Point4 multiplication
Point4 operator* (const Matrix4& m, const Point4& v) {
    // The "v" parameter is for "vector" (even though the type is Point4)
    Point4 ret;

    ret[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2] + m[3][0]*v[3];
    ret[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2] + m[3][1]*v[3];
    ret[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2] + m[3][2]*v[3];
    ret[3] = m[0][3]*v[0] + m[1][3]*v[1] + m[2][3]*v[2] + m[3][3]*v[3];

    return ret;
}

// Point * float scalar multiplication
// This works on Point2, Point3, and Point4
template <class T>
T operator* (float scalar, const T& vec) {
    T ret;
    for (int i = 0; i < vec.v.size(); i++) {
        ret.v[i] = vec.v[i] * scalar;
    }
    return ret;
}

template <class T>
T operator* (const T& vec, float scalar) {
    T ret;
    for (int i = 0; i < vec.v.size(); i++) {
        ret.v[i] = vec.v[i] * scalar;
    }
    return ret;
}

bool operator== (const Point2& a, const Point2& b) {
    return float_eq(a[0], b[0]) &&
           float_eq(a[1], b[1]);
}

bool operator== (const Point3& a, const Point3& b) {
    return float_eq(a[0], b[0]) &&
           float_eq(a[1], b[1]) &&
           float_eq(a[2], b[2]);
}

bool operator== (const Point4& a, const Point4& b) {
    return float_eq(a[0], b[0]) &&
           float_eq(a[1], b[1]) &&
           float_eq(a[2], b[2]) &&
           float_eq(a[3], b[3]);
}

// Note: The code on which this is based (PtahRenderer) defines only Point types (i.e., no types named "vector")

// cout a vector (adapted from https://stackoverflow.com/questions/10750057/how-to-print-out-the-contents-of-a-vector)
template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '[';

    //std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));  // This was the original line; but apparently, ostream_iterator is available only since C++17? But we're not using that standard for this project. Because reasons)
    for (auto element : v) {
        out << element << ", ";
    }

    out << "\b\b]";
  }
  return out;
}


// cout point/vec-like objects (i.e., my Point2, Point3, and Point4 classes)
// NOTE: this would be well-served as a template function, but I couldn't get it to work properly because of ambiguous overloads. Also, I hate templates, mainly because I don't REALLY understand tepmlating...
std::ostream& operator<<(std::ostream& out, const Point2& pointObj) {
    std::string sep = "   ";
    out << "[ ";
    int i;
    for (i = 0; i < pointObj.v.size() - 1; i++) {
        out << pointObj.v[i] << sep;
    }
    out << pointObj.v[i];
    out << " ]";
    return out;
}

std::ostream& operator<<(std::ostream& out, const Point3& pointObj) {
    std::string sep = "   ";
    out << "[ ";
    int i;
    for (i = 0; i < pointObj.v.size() - 1; i++) {
        out << pointObj.v[i] << sep;
    }
    out << pointObj.v[i];
    out << " ]";
    return out;
}

std::ostream& operator<<(std::ostream& out, const Point4& pointObj) {
    std::string sep = "   ";
    out << "[ ";
    int i;
    for (i = 0; i < pointObj.v.size() - 1; i++) {
        out << pointObj.v[i] << sep;
    }
    out << pointObj.v[i];
    out << " ]";
    return out;
}

// cout a Matrix4 
std::ostream& operator<< (std::ostream& out, const Matrix4& mat) {
    // Print the elements of mat as they would appear on paper
    out << mat.m[0][0] << "   " << mat.m[1][0] << "   " << mat.m[2][0] << "   " << mat.m[3][0] << std::endl;
    out << mat.m[0][1] << "   " << mat.m[1][1] << "   " << mat.m[2][1] << "   " << mat.m[3][1] << std::endl;
    out << mat.m[0][2] << "   " << mat.m[1][2] << "   " << mat.m[2][2] << "   " << mat.m[3][2] << std::endl;
    out << mat.m[0][3] << "   " << mat.m[1][3] << "   " << mat.m[2][3] << "   " << mat.m[3][3] << std::endl;
    return out;
}



