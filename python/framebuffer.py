from common import Color
from copy import deepcopy

# Note: in Python3, we don't have to subclass object (i.e., the "new-style
# class"), because it is the default super class
class Framebuffer:
    def __init__(self, w=0, h=0):
        ## "public"
        self.width = int(w)
        self.height = int(h)

        ## "private"
        # array of Color objects
        self._colorBuffer = [ Color() ] * self.width * self.height
        # array of floats (0.0 - 1.0)
        self._depthBuffer = [ 1.0 ] * self.width * self.height


    def getDepth(self, x: int, y: int) -> float:
        """ Get the depth from the depth buffer at a given x,y location
            
            :param x: x coordinate
            :type x: int
            :param y: y coordinate
            :type x: int
            :returns: float (the depth, 0.0 - 1.0, inclusive)
        """
        return self._depthBuffer[y * self.width + x]

    def setDepth(self, x: int, y: int, depth: float):
        """ Set the depth in the depth buffer at a given x,y location
            
            :param x: x coordinate
            :type x: int
            :param y: y coordinate
            :type x: int
            :param depth: the depth value [0.0 - 1.0], inclusive
            :type depth: float

            Note: There is no error checking or clamping here (in a more
            robust engine, clamping might be a good idea
        """
        self._depthBuffer[y * self.width + x] = depth
        
    def setColor(self, x: int, y: int, color: Color):
        self._colorBuffer[y * self.width + x] = deepcopy(color)

    def clear(self, color: Color = Color(0.0)):
        """ Clear the frame buffer

            By default, the color is black (0,0,0)
            But other colors can be specified
        """
        for i in range(0, self.width * self.height):
            self._colorBuffer[i] = color
            self._depthBuffer[i] = 1.0


    def write(self, path: str) -> bool:
        # 18 byte header
        header = [0] * 18

        header[2]  = 2
        header[8]  = 8
        header[9]  = 0
        header[10] = 0
        header[11] = 0
        header[12] = self.width % 256
        header[13] = self.width // 256
        header[14] = self.height % 256
        header[15] = self.height // 256
        header[16] = 24             # bits per pixel
        header[17] = 0              # image descriptor

        pixeldata = []

        for c in self._colorBuffer:
            pixeldata.extend( [int( max(min(c.b*255, 255), 0) ),
                               int( max(min(c.g*255, 255), 0) ),
                               int( max(min(c.r*255, 255), 0) )] )
        data = bytes(header + pixeldata)

        try:
            with open(path, "wb") as fd:
                fd.write(data)
        except Exception as e:
            print("Failed to write image file: {}".format(e))
            traceback.print_exc()
            return False

        return True


