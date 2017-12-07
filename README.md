This is the source to the native component of zOBBTransform, which is a part
of zMayaTools.

This source is relatively large, and I'd like to eliminate this dependency at
some point, so I'm putting this in a separate repository to avoid bloating zMayaTools.
The built Python module is 23k, so just that will be committed to zMayaTools.

Building this is annoying:

- Python modules only build with VC2013, not 2015.
- Install Python 2.7 for Windows.
- Run SET VS90COMNTOOLS=%VS120COMNTOOLS% to work around "Unable to find vcvarsall.bat".
- From ../obb_transform_src, run:
  "C:\Program Files\Python27\python.exe" setup.py build --build-lib ..\obb_transform

