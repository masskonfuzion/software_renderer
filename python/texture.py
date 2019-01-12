import os
from common import Color
import sdl, sdl2.ext

class Texture(object):
    # TODO -- remove "object" subclass? It's the default in Python3?
    def __init__(self, name):
        super().__init__()
        self.path = os.path.normpath("{}{}{}{}".format("Resources", os.sep, name, ".tga"))
        self.w = 0
        self.h = 0
        self._colors = []

        self._initialize()

    def __getattribute__(self, index):
        """ Return a indexed attribute of the _colors array.

            The item at index is, itself, an array.
            So, in typical usage, users will refer to, e.g., texture[x][y],
            where texture[x] returns an array (the "x column"), which is y
            elements long
        """
        # Note: apparently, __getattribute__ is available only in "modern"
        # Python; in older versions, __getattr__ would be required.
        return self._colors[index]

    def _load_image(self, path):
        image_surface = sdl2.ext.load_image(path)
        return image_surface

    def _initialize(self):
        """ Do some initializtion work
        """
        image_surface = self._load_image(self.path)

        # Get width and height 
        self.w = image_surface.clip_rect.w
        self.h = image_surface.clip_rect.h

        # Get pixels
        pixels = sdl2.ext.pixels2d(image_surface)

        # Initialize self._colors array (2D array, for simplicity)
        self._colors = self._get_colors_from_pixels(pixels)
        
        
    def _get_colors_from_pixels(self, pixels):
        """ Create an array of Color objects from an sdl2 pixels array

            NOTE: this function assumes 32 bits per pixel
        """
        RED_MASK   = 0x00ff0000
        GREEN_MASK = 0x0000ff00
        BLUE_MASK  = 0x000000ff

        colors = []

        # Ye olde double "for" loop (there's no way to go faster than this)
        # Note: pixels are stored as x,y (i.e. columns, then rows)
        for x in pixels:
            # Append an empty list to the colors array
            colors.append([])

            for y in pixels[x]:
                # note: we have to shift because the bitmask leaves the color
                # component where it was (i.e., the RED component is held in
                # bits 16..23 (0-indexed)
                r = (pixels[x][y] & RED_MASK) >> 16
                g = (pixels[x][y] & GREEN_MASK) >> 8
                b = (pixels[x][y] & BLUE_MASK)

                colors[x].append( Color(r / 255.0, g / 255.0, b / 255.0) )

        return colors
