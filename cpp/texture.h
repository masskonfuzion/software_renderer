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
        width = 0;
        height = 0;

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

        uint8_t IDLength = uint8_t(header[0]);
        width = uint8_t(header[13])*256 + uint8_t(header[12]);
        height = uint8_t(header[15])*256 + uint8_t(header[14]);
        uint8_t pixelDepth = uint8_t(header[16]);

        uint32_t lengthImage = pixelDepth * width * height / 8;

        fd.seekg(18 + IDLength);   // TODO delete?  We shouldn't need to seek because we've already read in the header. The file ptr shold be where we need it
        // Read entire file into memory at once into a dynamically allocated
        // memory block. We could also read byte-by-byte (which we would do,
        // if we were concerned about large files choking RAM, but here we're
        // not concerned about that).
        char* content = new char[lengthImage];
        fd.read(content, lengthImage);

        // Apparently our pointer type for binary file input must be char* (not uint8_t* a.k.a. unsigned char*)
        // So we're going to initialize another pointer, which treats our buffer as uint8_t
        // We're initializing our ptr so that the ptr itself is non-const (i.e. we can change what we're pointing at), but the data being pointed to is treated as const (i.e. the ptr is promising not to change the data)
        uint8_t* const content_u = reinterpret_cast<uint8_t* const> (content);

        if (pixelDepth == 8) {
            for (int i = 0; i < lengthImage; i++) {
                float byt = float(content_u[i]);
                pixels.push_back( Color(byt/255.0, byt/255.0, byt/255.0)  );
            }
        }
        else {
            int pixelBytes = pixelDepth / 8;

            for (int i = 0; i < width * height; i++) {
                float r = float(content_u[pixelBytes * i + 2]) / 255.0;
                float g = float(content_u[pixelBytes * i + 1]) / 255.0;
                float b = float(content_u[pixelBytes * i    ]) / 255.0;
                pixels.push_back(Color(r, g, b));
            }
        }
        delete[] content;   // Free memory used by texture image file
        // note that we delete content because it is actually the ptr that owns the data; however, now content AND content_u are invalid. Trying to dereference either at this point would cause a segfault
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
