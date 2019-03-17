#include <iostream>
#include "common.h"
#include "mesh.h"
#include "texture.h"
#include "rendering.h"
#include "clipping.h"
#include "framebuffer.h"

////TODO delete
//#include <string>
//#include <vector>

int main(int argc, char* argv[]) {
    std::cout << "Loading data..." << std::endl;

    //// TODO <delete>
    //std::vector<std::string> v;
    //std::string videogames = "nintendo sega sony microsoft";
    //string_split(videogames, ' ', v);
    //std::cout << v << std::endl;
    //// </delete>

    Mesh mesh("dragon");
    Texture texture("dragon");

    std::cout << "Preparing the scene..." << std::endl;

    float width = 800.0;
    float height = 600.0;

    Framebuffer framebuffer = Framebuffer(int(width), int(height));

    // Compute model-view-projection matrix (angle in radians);
    Matrix4 modelMatrix = Matrix4::scaleMatrix(0.5) * Matrix4::rotationMatrix(-1.7, Point3(0.0, 1.0, 0.0));

    Matrix4 viewMatrix = Matrix4::lookAtMatrix(Point3(-2.0, 0.5, 0.0), Point3(), Point3(0.0, 1.0, 0.0));

    Matrix4 projectionMatrix = Matrix4::perspectiveMatrix(70.0, width / height, 0.5, 10.0);

    // Create the MVP matrix (remember, we compose the transformations in
    // "reverse"
    Matrix4 mvpMatrix = projectionMatrix * viewMatrix * modelMatrix;

    Matrix4 lightMatrix = Matrix4::inverse( Matrix4::transpose( modelMatrix ) );

    std::cout << "Rendering..." << std::endl;

    // Fill the framebuffer with a constant background color
    //std::cout << "framebuffer size: (" << framebuffer.width << ", " << framebuffer.height << ")" << std::endl;
    framebuffer.clear( Color(0.1, 0.1, 0.1) );

    // For each face in the mesh...
    for (std::vector<Face>::iterator faceModelSpace = mesh.faces.begin(); faceModelSpace != mesh.faces.end(); faceModelSpace++) {
        // Transform the face into clip space
        Face faceClipSpace = worldSpaceToClipSpace(*faceModelSpace, mvpMatrix, lightMatrix);

        // Perform clipping, potentially producing additional faces
        std::vector<Face> clippedFaces = clip(faceClipSpace);

        for (std::vector<Face>::iterator clippedFace = clippedFaces.begin(); clippedFace != clippedFaces.end(); clippedFace++) {
            // Apply perspective division
            Face faceNDSpace = perspectiveDivide(*clippedFace);

            // check if the face is inviaible and should be culled
            if (cullFace(faceNDSpace)) {
                continue;
            }

            // Transform the face into screen space
            Face faceScreenSpace = normalizedSpaceToScreenSpace(faceNDSpace, width, height);

            // Draw the face
            // Note: the draw() function assumes a single texture -- in a
            // for-reals renderer, I'd probably associate the texture with a
            // face (the Face might have a pointer to a Material)
            draw(faceScreenSpace, framebuffer, texture);
        }
    }

    std::cout << "Saving the image..." << std::endl;
    bool success = framebuffer.write("result.tga");

    //if (success) {
    //    std::cout << "Success!";
    //}
    //else {
    //    std::cout << "Error while saving result"
    //}
    success ? std::cout << "Success!" : std::cout << "Error while saving result";

    std::cout << std::endl;

    return 0;
}
