#pragma once
#include <algorithm>
#include <fstream>
#include <vector>
#include "common.h"


struct Framebuffer {
public:
    Framebuffer(int w_in = 0, int h_in = 0) : width(w_in), height(h_in) {
        // Initialize the frame buffer
        colorBuffer_.clear();
        for (int i = 0; i < width * height; i++) {
            colorBuffer_.push_back(Color());
        }

        depthBuffer_.clear();
        for (int i = 0; i < width * height; i++) {
            depthBuffer_.push_back(1.0);
        }
    }

    // public member functions
    float getDepth(int x, int y) { return depthBuffer_[y * width + x]; }
    void setDepth(int x, int y, float depth) { depthBuffer_[y * width + x] = depth; }
    
    void setColor(int x, int y, Color color) { colorBuffer_[y * width + x] = color; }
    // ^^ setColor should copy the color object (C++ params are by-value (copy), by default)

    void clear(const Color& color) {
        for (int i = 0; i < width * height; i++) {
            // the 'brute-force'/'ugly' way to reinitialize colors in the color buffer
            colorBuffer_[i].r = color.r;
            colorBuffer_[i].g = color.g;
            colorBuffer_[i].b = color.b;

            // reset the depth
            depthBuffer_[i] = 1.0;
        }
    }

    bool write(std::string path) {
        std::ofstream fd;
        fd.open(path, std::ios::binary);
        bool status = false;    // return status

        if (fd.is_open()) {
            char header[18];
            // Zero the file header in an old-school, but explicit, manner
            // This isn't strictly necessary (it should be necessary only to
            // say char header[18] = {0x00}; but the behavior of an initialization
            // like that is less clear/well-defined in C++ than it is in other languages
            // e.g., int[5] iarray = {1}; actually yields [1, 0, 0, 0, 0]
            for (int i = 0; i < 18; i++) {
                header[i] = 0x00;
            }

            // Now, actually initialize the header
            header[2] = 0x02;
            header[8] = 0x08;
            header[12] = int(width) % 256;
            header[13] = int(width) / 256;
            header[14] = int(height) % 256;
            header[15] = int(height) / 256;
            header[16] = 0x18;              // bits per pixel (in hex. e.g. int 24 == char 0x18)
            header[17] = 0x00;              // image descriptor
            // write the header
            // (the size is sizeof(header), which works, because header is a
            // raw array)
            fd.write(header, sizeof(header));

            // Now, create the pixel data
            // Note: I wanted to do this as a std::array, but I got errors re: non-const expression,
            // because you can't create an array to a size that is not known at compile time (even if
            // the size you eventually want to use is of type "const xyz")
            //std::vector<uint8_t> pixeldata;   // TODO delete this line
            char pixeldata[width * height * 3];

            // TODO delete this for loop?
            //for (std::vector<Color>::iterator cbItem = colorBuffer_.begin(); cbItem != colorBuffer_.end(); ++cbItem) {
            //    // For each Color object, append 3 items to pixeldata (the float values for r, g, and b)
            //    pixeldata.push_back( std::max(std::min(int(cbItem.b*255), 255), 0) );
            //    pixeldata.push_back( std::max(std::min(int(cbItem.g*255), 255), 0) );
            //    pixeldata.push_back( std::max(std::min(int(cbItem.r*255), 255), 0) );
            //}

            for (int i = 0; i < colorBuffer_.size(); i++) {
                Color& cbItem = colorBuffer_[i];
                // For each Color object, write 3 items to pixeldata (the float values for r, g, and b)
                // Remember that the size of pixeldata is the size of the colorBuffer * 3
                // i.e., we store 3 separate ints (r, g, and b) for each Color
                pixeldata[i]   = std::max(std::min(int(cbItem.b*255), 255), 0);
                pixeldata[i+1] = std::max(std::min(int(cbItem.g*255), 255), 0);
                pixeldata[i+2] = std::max(std::min(int(cbItem.r*255), 255), 0);
            }
            fd.write(pixeldata, sizeof(pixeldata));

            fd.close();
            status = true;
        }
        else {
            status = false;
        }
        return status;
    }


    // public data members
    int width;
    int height;

private:
    // private data members
    // Note: these could also be std::array ptr (a ptr to a "new" (dynamically allocated) fixed-size array.
    // But the vector is easier to work with
    std::vector<Color> colorBuffer_;
    std::vector<float> depthBuffer_;
};
