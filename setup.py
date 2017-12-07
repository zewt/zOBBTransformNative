from distutils.core import setup, Extension

setup(
    name='obb_transform',
    version='1',
    py_modules=[],
    ext_modules=[Extension('obb_transform_native', ['obb_transform_native.cpp'])],
)
