#pragma once
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "common.h"

struct Texture {
    // Initialize a texture from a tga file
    // Note: This is hand-coded. With a more robust renderer, we could
    // probably use an image loader 
    Texture(std::string name) {
        // TODO use the std::filesystem classes/methods to determine how to construct the path (similar to using the os module in Python
        std::string path = "Resources" + os_path_sep + name + ".tga";

        char header[18];

        std::ifstream fd;
        fd.open(path, std::ios::binary);
        try
        {
            // Read 18 bytes into the header array
            fd.read(header, 18);
        }
        catch(std::exception& e)
        {
            // In this block, we catch any exception thrown
            // Our code doesn't actually throw any exceptions... But maybe
            // other code (e.g. STL) does
            std::cout << "Couldn't load the tga: " << e.what() << std::endl;
        }

        bool useColorMap = header[0] != 0x00;
        uint8_t imageType = uint8_t(header[2]);

        if (useColorMap ||
            imageType == 0x00 ||
            imageType == 0x01 ||
            imageType == 0x03)
        {
            // throw exception here
            throw std::runtime_error(path + " is not a color tga");
        }
    }

    // Public member functions
    // Accessor
    //Get the Color at a given u,v location
    //:param pos: Position in Texture image to return Color from
    const Color& operator[] (std::array<float, 2> pos) const {
        float a = pos[0] * width;
        float b = pos[1] * height;
        return subscript(int(a), int(b));
    }


    // Public data members
    int width;
    int height;
    std::vector<Color> pixels;

private:
    //Get the Color at a given x,y location
    //x,y is computed from u,v (normalized) texture coordinates
    const Color& subscript(int a, int b) const {
        return pixels[ std::min(std::max(b, 0), height-1) * width + std::min(std::max(a,0), width-1) ];
    }
};
