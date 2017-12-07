#define EIGEN_MPL2_ONLY

#include <Python.h>
#include <vector>

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

static PyObject *MatrixToPythonArray(const Eigen::MatrixXd &mat)
{
    PyObject *result = PyList_New(mat.rows());
    for(int y = 0; y < mat.rows(); ++y)
    {
        PyObject *row = PyTuple_New(mat.cols());
        for(int x = 0; x < mat.cols(); ++x)
            PyTuple_SET_ITEM(row, x, PyFloat_FromDouble(mat(y, x)));
        
        PyList_SET_ITEM(result, y, row);
    }
    return result;
}

static PyObject *VectorToPythonArray(const Eigen::VectorXd &vec)
{
    PyObject *result = PyList_New(vec.rows());
    for (int i = 0; i < vec.rows(); ++i)
        PyList_SET_ITEM(result, i, PyFloat_FromDouble(vec(i)));
    return result;
}

struct result
{
    void compute(const vector<Vector3d> &points);

    Vector3d right, up, forward;

    Vector3d m_rot[3];
    Vector3d m_pos;
    Vector3d m_ext;
};

void result::compute(const vector<Vector3d> &points)
{
    // loop over the points to find the mean point location
    Vector3d mu(0,0,0);
    for(int i=0; i< points.size(); i++)
        mu += points[i];
    mu /= double(points.size());

    // loop over the points again to build the covariance matrix.  Note that we only have
    // to build terms for the upper trianglular portion since the matrix is symmetric
    double cxx=0.0, cxy=0.0, cxz=0.0, cyy=0.0, cyz=0.0, czz=0.0;
    for(auto v: points)
    {
        v -= mu;
        cxx += v(0)*v(0); // - mu(0)*mu(0);
        cxy += v(0)*v(1); // - mu(0)*mu(0);
        cxz += v(0)*v(2); // - mu(0)*mu(0);
        cyy += v(1)*v(1); // - mu(1)*mu(1);
        cyz += v(1)*v(2); // - mu(1)*mu(1);
        czz += v(2)*v(2); // - mu(2)*mu(2);
    }

    // now build the covariance matrix
    Eigen::Matrix3d C;
    C(0,0) = cxx; C(0,1) = cxy; C(0,2) = cxz;
    C(1,0) = cxy; C(1,1) = cyy; C(1,2) = cyz;
    C(2,0) = cxz; C(2,1) = cyz; C(2,2) = czz;

    // extract the eigenvalues and eigenvectors from C
    SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(C);
    auto eigval = eigen.eigenvalues();
    auto eigvec = eigen.eigenvectors();

    // find the right, up and forward vectors from the eigenvectors
    forward = Vector3d(eigvec(0,0), eigvec(1,0), eigvec(2,0));
    up      = Vector3d(eigvec(0,1), eigvec(1,1), eigvec(2,1));
    right   = Vector3d(eigvec(0,2), eigvec(1,2), eigvec(2,2));
    right = right.normalized(); up = up.normalized(); forward = forward.normalized();

    // now build the bounding box extents in the rotated frame
    Vector3d minim(1e10, 1e10, 1e10), maxim(-1e10, -1e10, -1e10);
    for(const auto &point: points)
    {
        Vector3d p(forward.dot(point), up.dot(point), right.dot(point));
        minim = minim.cwiseMin(p);
        maxim = maxim.cwiseMax(p);
    }

    m_ext = (maxim-minim) * 0.5;

    // set the center of the OBB to the average of the minimum and maximum, and 
    // the extents to half of the difference between the minimum and maximum
    Vector3d center = (maxim+minim) * 0.5;

    // set the rotation matrix using the eigvenvectors
    m_rot[0](0)=forward(0); m_rot[0](1)=up(0); m_rot[0](2)=right(0);
    m_rot[1](0)=forward(1); m_rot[1](1)=up(1); m_rot[1](2)=right(1);
    m_rot[2](0)=forward(2); m_rot[2](1)=up(2); m_rot[2](2)=right(2);

    m_pos = Vector3d( m_rot[0].dot(center), m_rot[1].dot(center), m_rot[2].dot(center) );
}

static PyObject *obb_transform(PyObject *self, PyObject *args)
{
    result r;
    vector<Vector3d> points;
    int rank = 10;

    PyObject *matrix_arg = NULL;

    if(!PyArg_ParseTuple(args, "O|i", &matrix_arg, &rank))
        goto error;
    if(rank < 1)
        rank = 1;

    PyObject *ret = NULL;

    PyObject *matrix = PySequence_Fast(matrix_arg, "expected a sequence");
    if(PyErr_Occurred())
        goto error;

    int rows = (int) PySequence_Size(matrix_arg);
    if(rows == 0)
    {
        PyErr_SetString(PyExc_TypeError, "Empty matrix");
        goto error;
    }

    points.resize(rows);
    for(int y = 0; y < rows; ++y)
    {
        PyObject *row = PySequence_Fast_GET_ITEM(matrix, y);
        PyObject *row_data = PySequence_Fast(row, "expected a sequence");
        if(PyErr_Occurred())
            goto error;

        int cols = (int) PySequence_Size(row_data);

        if(cols != 3)
        {
            PyErr_SetString(PyExc_TypeError, "Inconsistent column count");
            Py_XDECREF(row_data);
            goto error;
        }

        for(int x = 0; x < 3; ++x)
        {
            PyObject *value = PySequence_Fast_GET_ITEM(row_data, x);
            points[y](x) = PyFloat_AsDouble(value);
        }

        Py_XDECREF(row_data);

        if(PyErr_Occurred())
            goto error;
    }

    r.compute(points);

    PyObject *forward = VectorToPythonArray(r.right); // +x
    PyObject *up = VectorToPythonArray(r.up); // +y
    PyObject *right = VectorToPythonArray(r.forward); // +z
    PyObject *center = VectorToPythonArray(r.m_pos);
    PyObject *extent = VectorToPythonArray(r.m_ext);

    ret = PyTuple_Pack(5, forward, up, right, center, extent);

error:
    Py_XDECREF(matrix);
    
    return (PyObject *)ret;
}

static PyMethodDef obb_transform_native_methods[] = {
    {"obb_transform", (PyCFunction)obb_transform, METH_VARARGS, ""},
    {NULL, NULL},
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "obb_transform_native",
    NULL,
    0,
    obb_transform_native_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};
#endif

#if PY_MAJOR_VERSION >= 3
#define INIT_ERROR return NULL
#else
#define INIT_ERROR return
#endif

#if PY_MAJOR_VERSION >= 3
PyObject *PyInit__obb_transform_native()
{
    PyObject *module = PyModule_Create(&moduledef);
#else
void initobb_transform_native()
{
    PyObject *module = Py_InitModule("obb_transform_native", obb_transform_native_methods);
#endif
    if (module == NULL) {
        INIT_ERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

/* 
 *  Copyright (c) 2010 Daisuke Okanohara
 *  Copyright (c) 2016 Glenn Maynard
 * 
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 * 
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

