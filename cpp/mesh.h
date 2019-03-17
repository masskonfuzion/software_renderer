#pragma once

//#include <filesystem>   // added as of C++17 -- but do I need to compile with '-std=c++17'?
#include <exception>
#include <fstream>
#include <array>
#include <string>
#include <vector>
#include "common.h"

struct Vertex {
    // Default constructor
    Vertex(){}
    // Overloaded default constructor
    Vertex(const Point4& p_in, const Point4& n_in, const Point2& t_in) : p(p_in), n(n_in), t(t_in) {}
    //// Copy constructor
    Vertex(const Vertex& src) {
        p = src.p;
        n = src.n;
        t = src.t;
    }


    // data members
    Point4 p;
    Point4 n;
    Point2 t;
};

// A free-floating equality operator
bool operator== (const Vertex& a, const Vertex& b) {
    return a.p == b.p &&
           a.n == b.n &&
           a.t == b.t;
}


struct Face {
    // Default constructor
    Face(){}
    // Overloaded default constructor (hopefully uses copy constructors to copy the Vertexes)
    // TODO const correctness for Face constructor? (we'll need to make sure we're deep-copying the passed-in parameters)
    Face(const Vertex& v0_in, const Vertex& v1_in, const Vertex& v2_in) : v0(v0_in), v1(v1_in), v2(v2_in) {}
    // A second overloaded default constructor
    // (to avoid manually replicating unchanged normals and texture coordinates
    // when updating positions)
    // TODO const correctness for Face constructor?
    Face(const Point4& p0, const Point4& p1, const Point4& p2, const Face& face) {
        v0 = Vertex(p0, face.v0.n, face.v0.t);
        v1 = Vertex(p1, face.v1.n, face.v1.t);
        v2 = Vertex(p2, face.v2.n, face.v2.t);
    }

    // Copy constructor
    Face(const Face& src) {
        v0 = src.v0;
        v1 = src.v1;
        v2 = src.v2;
    }

    Vertex v0;
    Vertex v1;
    Vertex v2;
};


struct Mesh {
    // Default constructor
    Mesh(){}
    // Overloaded default constructor
    Mesh(std::string name) {
        struct FaceIndices {
            // Default constructor
            FaceIndices() : p(0), t(0), n(0) {}
            // TODO -- make a copy constructor? We might not NEED to--C++ gives us a copy constructor that may work for this type
            // Overloaded default constructor
            FaceIndices(int p_in, int t_in, int n_in) : p(p_in), t(t_in), n(n_in) {}
            int p;
            int n;
            int t;
        };

        // TODO use the std::filesystem classes/methods to determine how to construct the path (similar to using the os module in Python
        std::string path = "Resources" + os_path_sep + name + ".obj";
        //std::cout << "DEBUG: Created a Mesh object with path: " << path << std::endl;

        std::vector<Point4> positions;
        std::vector<Point4> normals;
        std::vector<Point2> uvs;
        std::vector< std::array<FaceIndices, 3> > indices;


        // TODO try/catch/throw or whatever
        try
        {
            std::ifstream fd;
            std::string line_in;
            fd.open(path);
            if(fd.is_open()) {  // Note: this check feels unnecessary.. But it's here for good measure
                std::string lineprefix;
                std::vector<std::string> components; // i.e., components of each line item
                while(std::getline(fd, line_in)) {
                    string_split(line_in, ' ', components);
                    //std::cout << line_in << std::endl;  // debug print the line as read from file
                    std::cout << components << std::endl;   // debug print the line as split into components

                    lineprefix = components[0];

                    if (lineprefix == "v") {
                        // Position (Point4: x,y,z,1)
                        positions.push_back( Point4(std::stof(components[1]),
                                                     std::stof(components[2]),
                                                     std::stof(components[3])) );
                    }
                    else if (lineprefix == "vt") {
                        // UV texture coordinates
                        uvs.push_back( Point2(std::stof(components[1]),
                                              std::stof(components[2])) );
                    }
                    else if (lineprefix == "vn") {
                        // Normal coordinates (Point4: x,y,z,0)
                        normals.push_back( Point4(std::stof(components[1]),
                                                  std::stof(components[2]),
                                                  std::stof(components[3])) );
                    }
                    else if (lineprefix == "f") {
                        // Rearrange a string, like 
                        //   "f  v1/vt1/vn1   v2/vt2/vn2   v3/vt3/vn3"
                        // into list of lists of strings, like 
                        //   [ [v1, vt1, vn1], [v2, vt2, vn2], [v3, vt3, vn3] ]
                        // Example:
                        //   "f 7/1/1 8/2/2 9/3/3"
                        // turns into
                        //   [ ["7", "1", "1"], ["8", "2", "2"], ["9", "3", "3"] ]

                        std::vector<std::string> tmp_indices;   // strings (vector required by string_split())
                        std::array<std::array<int, 3>, 3> ids_by_vtx;   // An array of arrays

                        // Iterate of the line components (i.e. the substr's that were split on ' '),
                        // ignoring the lineprefix
                        for (int i = 1; i < components.size(); i++) {
                            string_split(components[i], '/', tmp_indices);

                            // Convert indices from string to int, store in ids_by_vtx
                            for (int j = 0; j < tmp_indices.size(); j++) {
                               ids_by_vtx[i-1][j] = std::stoi(tmp_indices[j]);
                            }
                        }

                        std::array<FaceIndices, 3> fi;  // ids_by_vtx.size() below == 3; that's why this works
                        for (int i = 0; i < ids_by_vtx.size(); i++) {
                            // FaceIndices constructor takes parameters: p (Position index), t (Tex/UV) index, n (Normal) index
                            fi[i] = FaceIndices(ids_by_vtx[i][0]-1, ids_by_vtx[i][1]-1, ids_by_vtx[i][2]-1);
                        }
                        indices.push_back(fi); // I think this should work, because push_back copies (or moves) the items being passed in
                    }
                }

                // iterate over the indices vector -- construct faces
                // TODO replace counter-based loop with iterator
                for (int i = 0; i < indices.size(); i++) {
                    // TODO 2019-03-02 - continue from here! Fill in the indices
                    faces.push_back( Face ( Vertex( positions[ indices[i][0].p ], normals[ indices[i][0].n ], uvs[ indices[i][0].t ]),
                                            Vertex( positions[ indices[i][1].p ], normals[ indices[i][1].n ], uvs[ indices[i][1].t ]),
                                            Vertex( positions[ indices[i][2].p ], normals[ indices[i][2].n ], uvs[ indices[i][2].t ])) );
                }

                //for (std::vector<std::array<FaceIndices, 3>>::iterator idx_obj = indices.begin(); idx_obj != indices.end(); idx_obj++) {
                //    faces.push_back( Face ( Vertex
                //}

                fd.close();
            }
            else {
                // Do something, maybe? We reach this block if the file is is not open (i.e., fd.is_open() returned false)
            }
        }
        catch(std::exception& e)
        {
            // In this block, we catch any exception thrown
            // Our code doesn't actually throw any exceptions... But maybe
            // other code (e.g. STL) does
            std::cout << "Couldn't load the mesh: " << e.what() << std::endl;
        }
    }

    // Copy constructor (will we be copying meshes, though?)
    Mesh(const Mesh& src) {
        faces.clear();
        for (std::vector<Face>::const_iterator face = src.faces.begin(); face != src.faces.end(); face++) {
            faces.push_back(*face); // Copy src's faces vector, one face at a time
        }
    }

    // Data members
    std::vector<Face> faces;
};
