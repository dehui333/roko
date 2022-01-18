#include <Python.h>
#include <cstdio>
#include <iostream>

#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"

#include "generate.h"

// Module method definitions
static PyObject* generate_features_cpp(PyObject *self, PyObject *args) {
    srand(time(NULL));

    char *filename, *ref, *region;
    PyObject* dict = NULL;
    int use_dict = -1;

    if (!PyArg_ParseTuple(args, "sssOi", &filename, &ref, &region, &dict, &use_dict)) return NULL;
    auto map = convert_py_labels_dict(dict);	
    auto result = generate_features(filename, ref, region, map, use_dict);

    PyObject* return_tuple = PyTuple_New(3);
    PyObject* pos_list = PyList_New(result->positions.size());
    PyObject* X_list = PyList_New(result->X.size());
    PyObject* Y_list = PyList_New(result->Y.size());

    for (int i = 0, size=result->positions.size(); i < size; i++) {
        auto& pos_element = result->positions[i];

        PyObject* inner_list = PyList_New(pos_element.size());
        for (int j = 0, s = pos_element.size(); j < s; j++) {
            PyObject* pos_tuple = PyTuple_New(2);
            PyTuple_SetItem(pos_tuple, 0, PyLong_FromLong(pos_element[j].first));
            PyTuple_SetItem(pos_tuple, 1, PyLong_FromLong(pos_element[j].second));

            PyList_SetItem(inner_list, j, pos_tuple);
        }
        PyList_SetItem(pos_list, i, inner_list);

        PyList_SetItem(X_list, i, result->X[i]);
        PyList_SetItem(Y_list, i, result->Y[i]);
    }
 
    PyTuple_SetItem(return_tuple, 0, pos_list);
    PyTuple_SetItem(return_tuple, 1, X_list);
    PyTuple_SetItem(return_tuple, 2, Y_list);

    return return_tuple;
}

static PyMethodDef gen_methods[] = {
        {
                "generate_features", generate_features_cpp, METH_VARARGS,
                "Generate features for polisher."
        },
        {NULL, NULL, 0, NULL}
};


static struct PyModuleDef gen_definition = {
        PyModuleDef_HEAD_INIT,
        "gen",
        "Feature generation.",
        -1,
        gen_methods
};


PyMODINIT_FUNC PyInit_gen(void) {
    Py_Initialize();
    import_array();
    return PyModule_Create(&gen_definition);
}
