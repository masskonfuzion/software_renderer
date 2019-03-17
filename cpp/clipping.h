#pragma once
#include <cmath>
#include <vector>
#include "mesh.h"

void clipEdge(const Vertex& v0, const Vertex& v1, std::vector<Vertex>& vertices) {
    // TODO determine whether or not we need to write a copy constructor for Vertex (I think the default _should_ work, but I'm not sure
    Vertex v0New = v0;
    Vertex v1New = v1;

    // Testing with respect to near plane
    bool v0Inside = v0.p.w > 0.0 && v0.p.z > -v0.p.w;
    bool v1Inside = v1.p.w > 0.0 && v1.p.z > -v1.p.w;

    if (v0Inside && v1Inside) {
        // Great! Nothing to do :-D
    }
    else if (v0Inside || v1Inside) {
        // Compute interpolation coefficients
        float d0 = v0.p.z + v0.p.w;
        float d1 = v1.p.z + v1.p.w;
        float factor = 1.0 / (d1 - d0);

        // New vertex with interpolated coefficients
        Vertex newVertex( factor * (d1 * v0.p - d0 * v1.p),
                          factor * (d1 * v0.n - d0 * v1.n),
                          factor * (d1 * v0.t - d0 * v1.t) );

        if (v0Inside) {
            v1New = newVertex;
        }
        else {
            v0New = newVertex;
        }
    }
    else {
        // Both are outside, on the same side. Remove the edge
        return;
    }


    // Add the first vertex if not already added
    int num_vertices = vertices.size();
    if (num_vertices == 0 || !(vertices[num_vertices-1] == v0New)) {
        // "not (v1 == v2)", instead of "v1 != v2", because we never
        // implemented a "!=" function.. I know... lazy :-)
        vertices.push_back(v0New);
    }

    // Add the second vertex
    vertices.push_back(v1New);
}


std::vector<Face> clip(const Face& faceClipSpace) {
    std::vector<Face> triangles;

    // All vertices are behind the camera
    if (faceClipSpace.v0.p.w <= 0.0 &&
        faceClipSpace.v1.p.w <= 0.0 &&
        faceClipSpace.v2.p.w <= 0.0 ) {

        // empty vector
        // TODO -- make sure we use move semantics when returning objects from functions?
        return triangles;
    }

    // All vertices are in front of the camera and inside the frustum
    if (faceClipSpace.v0.p.w > 0.0 && faceClipSpace.v1.p.w > 0.0 && faceClipSpace.v2.p.w > 0.0 &&
        std::abs(faceClipSpace.v0.p.z) < faceClipSpace.v0.p.w &&
        std::abs(faceClipSpace.v1.p.z) < faceClipSpace.v1.p.w &&
        std::abs(faceClipSpace.v2.p.z) < faceClipSpace.v2.p.w) {

        triangles.push_back(faceClipSpace);     // A list (vector) with a single face in it
    }
    else {
        // Clip each edge, accumulating vertices that we add or keep in an array
        std::vector<Vertex> vertices;

        clipEdge(faceClipSpace.v0, faceClipSpace.v1, vertices);
        clipEdge(faceClipSpace.v1, faceClipSpace.v2, vertices);
        clipEdge(faceClipSpace.v2, faceClipSpace.v0, vertices);

        // If not enough vertices to create a triangular face
        if (vertices.size() < 3) {
            triangles.clear();
            return triangles;
        }

        // We potentially have a duplicated at the end, that we can remove
        if (vertices[vertices.size()-1] == vertices[0]) {
            vertices.pop_back();    // Remove from the end
        }

        // Generate a fan of triangles, all sharing the first vertex
        for (int i = 1; i < vertices.size()-1; i++) {
            triangles.push_back(Face(vertices[0], vertices[i], vertices[i+1]));
        }
    }

    return triangles;
}


bool cullFace(const Face& faceNormalizedSpace) {
    float d = (faceNormalizedSpace.v1.p.x - faceNormalizedSpace.v0.p.x) *
              (faceNormalizedSpace.v2.p.y - faceNormalizedSpace.v0.p.y) -
              (faceNormalizedSpace.v1.p.y - faceNormalizedSpace.v0.p.y) *
              (faceNormalizedSpace.v2.p.x - faceNormalizedSpace.v0.p.x);
    return d < 0.0;
}

