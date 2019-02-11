# MKF Renderer

Some important specs and what not:

Python version:
* Python 3

Dependencies:
* numpy for math

This code is based on https://github.com/kosua20/PtahRenderer.git (see also the blog post at http://blog.simonrodriguez.fr/articles/18-02-2017_writing_a_small_software_renderer.html).
* The code is a Python port of the original code, which done in Swift
* The Resources dir is a direct copy of the "demo" branch of the PtahRenderer repository

This renderer writes a single frame to a file.  And, it is _unconscionably_ slow, like, taking 90 sec to render the dragon image and save it to disk, on a not-too-shabby Intel Core i7 at 2.4 GHz (quad-core, hyperthreading enabled).  In future work, I'd like to do a real-time renderer.
