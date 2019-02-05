from common import Color
import struct
import os
import traceback

class Texture:
    def __init__(self, name):
        """ Initialize a texture from a tga file
            
            Note: This is hand-coded. With a more robust renderer, we could
            probably use an image loader 
        """
        self._width : int = 0
        self._height : int = 0
        self._pixels : [Color] = []

        path = "{}{}{}{}".format('Resources', os.path.sep, name, '.tga')

        with open(path, 'rb') as fd:
            try:
                # Read file header
                # TODO worry about sign? (we should be reading uint8)
                # Note: header is of type "bytes"
                header = fd.read(18)
            except Exception as e:
                print("Couldn't load the tga {}: {}".format(path, e))
                traceback.print_exc()
                return

            useColorMap = header[0] != 0
            imageType = int(header[2])
            if useColorMap or imageType in (0, 1, 3):
                raise Exception("{} is not a color tga".format(path))
                return

            IDLength = header[0]
            self._width = int(header[13])*256 + int(header[12])
            self._height = int(header[15])*256 + int(header[14])
            pixelDepth = header[16]

            lengthImage = int(pixelDepth) * self._width * self._height // 8

            fd.seek(18 + int(IDLength))
            content = fd.read(lengthImage)

            if pixelDepth == 8:
                self._pixels = [ Color(float(byt)/255.0, float(byt)/255.0, float(byt)/255.0) for byt in content ]
            else:
                pixelBytes = pixelDepth // 8

                for i in range(0, (self._width * self._height)):
                    r = float(content[pixelBytes * i + 2]) / 255.0
                    g = float(content[pixelBytes * i + 1]) / 255.0
                    b = float(content[pixelBytes * i    ]) / 255.0
                    self._pixels.append(Color(r, g, b))


    ##def _subscript(self, a: int, b: int) -> Color:
    ##    return self._pixels

    ##def __getitem__(self, pos):
    ##    """ Get the Color at a given x,y location

    ##        :param pos: Position in Texture image to return Color from
    ##        :type pos: 2-tuple: float
    ##    """
    ##    a = pos[0] * float(self._width)
    ##    b = pos[1] * float(self._height)

    ##    return self._pixels[ min(max(b,0), height-1) * width + \
    ##                         min(max(a,0), width-1) ]


    def __getitem__(self, pos: tuple) -> Color:
        """ Get the Color at a given u,v location

            :param pos: Position in Texture image to return Color from
            :type pos: 2-tuple: float (values from 0.0 - 1.0)
        """
        a = pos[0] * self._width
        b = pos[1] * self._height
        return self._subscript( int(pos[0]), int(pos[1]))

    def _subscript(self, a: int, b: int) -> Color:
        """ Get the Color at a given x,y location

            x,y is computed from u,v (normalized) texture coordinates

            :param a: 
        """
        return self._pixels[ min(max(b, 0), self._height-1) * self._width + min(max(a, 0), self._width-1) ]
