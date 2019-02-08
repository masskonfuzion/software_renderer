#!/usr/bin/env python

from common import *
from mesh import *
from texture import *
from rendering import *
from clipping import *
from framebuffer import *



def main():
    print("Loading data...")

    mesh = Mesh("dragon")
    texture = Texture("dragon")
    
    print("Preparing the scene...")

    (width, height) = 800.0, 600.0
    framebuffer = Framebuffer(w=width, h=height)

    # Compute model-view-projection matrix (angle in radians)
    modelMatrix = Matrix4.scaleMatrix(0.5).dot(Matrix4.rotationMatrix(angle=-1.7, axis=Point3(0.0, 1.0, 0.0)))

    viewMatrix = Matrix4.lookAtMatrix(eye=Point3(-2.0, 0.5, 0.0), target=Point3(), up=Point3(0.0, 1.0, 0.0))

    projectionMatrix = Matrix4.perspectiveMatrix(fov=70.0, aspect=width/height, near=0.5, far=10.0)

    # Remember: the transformation is written PVM*p
    # Mathematically, transformations are applied in the following order:
    #   Model, View, Projection
    #mvpMatrix = projectionMatrix * viewMatrix * modelMatrix
    mvpMatrix = projectionMatrix.dot(viewMatrix.dot(modelMatrix))

    lightMatrix = numpy.linalg.inv(modelMatrix).transpose()


    print("Rendering...")

    # Fill the framebuffer with a constant background color
    framebuffer.clear(color=Color(0.1, 0.1, 0.1))

    # For each face in the mesh...
    for faceModelSpace in mesh.faces:
        
        # Transform the face into clip space
        faceClipSpace = worldSpaceToClipSpace(faceModelSpace, mvpMatrix, lightMatrix)

        # Perform clipping, potentially producing additional faces
        clippedFaces = clip(faceClipSpace)

        # For each of the clipped faces...
        for clippedFace in clippedFaces:

            # Apply perspective division
            faceNDSpace = perspectiveDivide(clippedFace)

            # Check if the face is invisible and ahould be culled
            if cullFace(faceNDSpace):
                continue

            # Transform the face into screen space
            faceScreenSpace = normalizedSpaceToScreenSpace(faceNDSpace, width, height)

            # Draw the face
            # Note: the draw() function assumes a single texture -- in a
            # for-reals renderer, I'd probably associate the texture with a
            # Face (the Face might contain or point to a Material)
            draw(faceScreenSpace, framebuffer, texture)
    
    print("Saving the image...")

    success = framebuffer.write("result.tga")

    print("Done!" if success else "Error while saving result")

if __name__ == "__main__":
    main()
