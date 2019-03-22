#pragma once
#include <cmath>
#include "common.h"
#include "framebuffer.h"
#include "mesh.h"
#include "texture.h"

Face worldSpaceToClipSpace(Face faceModelSpace, Matrix4 mvp, Matrix4 lightMatrix) {
    Vertex v0( Point4(mvp * faceModelSpace.v0.p),
               Point4(lightMatrix * faceModelSpace.v0.n),
               faceModelSpace.v0.t );
    Vertex v1( Point4(mvp * faceModelSpace.v1.p),
               Point4(lightMatrix * faceModelSpace.v1.n),
               faceModelSpace.v1.t );
    Vertex v2( Point4(mvp * faceModelSpace.v2.p),
               Point4(lightMatrix * faceModelSpace.v2.n),
               faceModelSpace.v2.t );

    return Face(v0, v1, v2);
}


// Transform model from clip space to normalized space
Face perspectiveDivide (const Face& faceClipSpace) {
    // TODO determine -- do we want to pass by reference here, and operate on the original face? Or do we want to copy?

    float w0_inv = 1.0 / faceClipSpace.v0.p.w;
    float w1_inv = 1.0 / faceClipSpace.v1.p.w;
    float w2_inv = 1.0 / faceClipSpace.v2.p.w;

    Point4 p0NormalizedSpace( faceClipSpace.v0.p.x * w0_inv,
                              faceClipSpace.v0.p.y * w0_inv,
                              faceClipSpace.v0.p.z * w0_inv,
                              faceClipSpace.v0.p.w );

    Point4 p1NormalizedSpace( faceClipSpace.v1.p.x * w1_inv,
                              faceClipSpace.v1.p.y * w1_inv,
                              faceClipSpace.v1.p.z * w1_inv,
                              faceClipSpace.v1.p.w );

    Point4 p2NormalizedSpace( faceClipSpace.v2.p.x * w2_inv,
                              faceClipSpace.v2.p.y * w2_inv,
                              faceClipSpace.v2.p.z * w2_inv,
                              faceClipSpace.v2.p.w );

    return Face(p0NormalizedSpace, p1NormalizedSpace, p2NormalizedSpace, faceClipSpace);
}


// Convert face from normalized space into screen space
Face normalizedSpaceToScreenSpace(const Face& faceNormalizedSpace, int width, int height) {
    
    Point4 p0ScreenSpace( floor(0.5 * width * (faceNormalizedSpace.v0.p.x + 1.0)),
                          floor(0.5 * height * (faceNormalizedSpace.v0.p.y + 1.0)),
                          faceNormalizedSpace.v0.p.z,
                          faceNormalizedSpace.v0.p.w );

    Point4 p1ScreenSpace( floor(0.5 * width * (faceNormalizedSpace.v1.p.x + 1.0)),
                          floor(0.5 * height * (faceNormalizedSpace.v1.p.y + 1.0)),
                          faceNormalizedSpace.v1.p.z,
                          faceNormalizedSpace.v1.p.w );

    Point4 p2ScreenSpace( floor(0.5 * width * (faceNormalizedSpace.v2.p.x + 1.0)),
                          floor(0.5 * height * (faceNormalizedSpace.v2.p.y + 1.0)),
                          faceNormalizedSpace.v2.p.z,
                          faceNormalizedSpace.v2.p.w );

    return Face(p0ScreenSpace, p1ScreenSpace, p2ScreenSpace, faceNormalizedSpace);
}


// Return a 2-tuple of Point2 objects, representing a a bounding box around screen space
std::array<Point2, 2> boundingBox(const Face& vs, int width, int height) {
    // We work only in the (x,y) plane.
    Point2 v0ScreenSpace(vs.v0.p.x, vs.v0.p.y);
    Point2 v1ScreenSpace(vs.v1.p.x, vs.v1.p.y);
    Point2 v2ScreenSpace (vs.v2.p.x, vs.v2.p.y);

    // Find minimal and maximal points.
    // Note -- these are pairwise min/max's
    Point2 mini = Point2::min( Point2::min(v0ScreenSpace, v1ScreenSpace), v2ScreenSpace );
    Point2 maxi = Point2::max( Point2::max(v0ScreenSpace, v1ScreenSpace), v2ScreenSpace );

    // Framebuffer bounds.
    Point2 lim( float(width) - 1.0, float(height) - 1.0 );

    // Clamp the bounding box against the framebuffer.
    Point2 finalMin = Point2::clamp( Point2::min(mini, maxi), Point2(0,0), lim );
    Point2 finalMax = Point2::clamp( Point2::max(mini, maxi), Point2(0,0), lim );
    return std::array<Point2, 2> {finalMin, finalMax};
}


Point3 barycentre(const Point2& p, const Point4& v0, const Point4& v1, const Point4& v2) {
	// v0 will be the origin.
	// ab and ac: basis vectors.
    Point4 ab = v1 - v0;
    Point4 ac = v2 - v0;

    // pa: vector with p coordinates in the barycentric frame.
    Point2 pa = Point2(v0.x, v0.y) - p;

    // magic
    Point3 uv1 = cross( Point3(ac.x, ab.x, pa.x), Point3(ac.y, ab.y, pa.y) );

    // Avoid division imprecision
    if (std::abs(uv1.z) < 1e-2) {
        return Point3(-1.0, 1.0, 1.0);
    }
	return (1.0/uv1.z)*Point3(uv1.z-(uv1.x+uv1.y), uv1.y, uv1.x);
}


// Update the frame buffer based on incoming Face data
void draw(const Face& fScreen, Framebuffer& framebuffer, const Texture& texture) {
    // compute clamped bounding box
    // boundingBox returns a tuple
    std::array<Point2, 2> bounds = boundingBox(fScreen, framebuffer.width, framebuffer.height);

    // TODO delete the next 2 lines?
    const Point2& mini = (const Point2&) bounds[0];
    const Point2& maxi = (const Point2&) bounds[1];

    for(int x = mini.x; x <= maxi.x; x++) {
        for (int y = mini.y; y <= maxi.y; y++) {
            // Compute the barycentric coordinates of the current pixel
            Point3 bary = barycentre( Point2(float(x), float(y)), fScreen.v0.p, fScreen.v1.p, fScreen.v2.p);

            // If one of them is negative, we are outside the triangle
            if (bary.x < 0.0 || bary.y < 0.0 || bary.z < 0.0) {
                continue;
            }

            // Interpolate depth at the current pixel
            float z = fScreen.v0.p.z * bary.x +
                      fScreen.v1.p.z * bary.y +
                      fScreen.v2.p.z * bary.z;

            // If the curent triangle pixel is closer than the last one drawn
            if (z < framebuffer.getDepth(x, y)) {
                // Compute perspective-correct interpolation coefficients
                Point3 persp( bary.x / fScreen.v0.p.w,
                              bary.y / fScreen.v1.p.w,
                              bary.z / fScreen.v2.p.w );

                // TODO write a *= operator for the vector types -- I feel like it will be a shade faster than creating a new variable on the stack
                persp = (1.0 / (persp.x + persp.y + persp.z)) * persp;

                // Perspective interpolation of texture coordinates and normal
                Point2 tex = persp.x * fScreen.v0.t + persp.y * fScreen.v1.t + persp.z * fScreen.v2.t;
                Point4 nor = persp.x * fScreen.v0.n + persp.y * fScreen.v1.n + persp.z * fScreen.v2.n;

                // Compute the color: use a diffuse factor with a strong light intensity
                float diffuse = 1.5 * std::max(float(0.0), dot(normalize(Point3(nor.x, nor.y, nor.z)), normalize(Point3(0.0, 1.0, 1.0))) );
                //Color color = diffuse * texture[ {tex.x, tex.y} ];
                Color color = diffuse * texture[ {tex.x, tex.y} ];

                // Set the framebuffer's color and depth at location (x,y)
                framebuffer.setColor(x, y, color);
                framebuffer.setDepth(x, y, z);
            }
        }
    }
}

