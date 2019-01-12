from common import Color

class Framebuffer(object):
    """ Depths are stored as normalized values (0.0 - 1.0)
        Colors are stored as objects with r, g, and b components
    """
    def __init__(self, width, height):
        super().__init__()
        self.w = width
        self.h = height
        self._depth_buffer = [1.0] * width * height
        self._color_buffer = [ Color(0.0, 0.0, 0.0) ] * width * height

    def getDepth(self, x, y):
        pass
        #return self._depth_buffer

    def setDepth(self, x, y, depth):
        pass

    def getColor(self, x, y):
        pass

    def setColor(self, x, y, color):
        pass

