#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <base/linalg.h>

namespace py = pybind11;

PYBIND11_MODULE(moopy, m)
{
    // ---  Linalg Namespace Bindings   ---
    py::module_ linalg_m = m.def_submodule("linalg", "Bindings for the Linalg C/C++ namespace.");

    // void square(int size, f64* x) -> x_cpy
    linalg_m.def("square", [](
        py::array_t<f64, py::array::c_style | py::array::forcecast> x_arr)
    {
        Linalg::square(static_cast<int>(x_arr.size()), x_arr.mutable_data());
        return x_arr;
    },
    py::arg("x"), "Perform element-wise squaring: y[i] := x[i] * x[i].");

    // f64 dot(int size, const f64* x, const f64* y)
    linalg_m.def("dot", [](
        py::array_t<const f64, py::array::c_style | py::array::forcecast> x_arr,
        py::array_t<const f64, py::array::c_style | py::array::forcecast> y_arr)
    {
        if (x_arr.size() != y_arr.size()) throw std::runtime_error("Vectors must have the same length for dot product.");
        return Linalg::dot(static_cast<int>(x_arr.size()), x_arr.data(), y_arr.data());
    },
    py::arg("x"), py::arg("y"), "Perform dot product: ret = transpose(x) * y.");

    // void matrix_vector(int size, char format, const double* matrix,
    //                    const double* vector, double* out) -> out
    linalg_m.def("matrix_vector", [](
        char format,
        py::array_t<const f64, py::array::c_style | py::array::forcecast> matrix_arr,
        py::array_t<const f64, py::array::c_style | py::array::forcecast> vector_arr)
    {
        if (matrix_arr.size() != (vector_arr.size() * vector_arr.size())) throw std::runtime_error("Matrix array size must be square of vector size.");

        py::array_t<f64> out_arr(vector_arr.size());
        Linalg::matrix_vector(
            vector_arr.size(), 
            format, 
            matrix_arr.data(), 
            vector_arr.data(), 
            out_arr.mutable_data()
        );
        return out_arr;
    },
    py::arg("format"), py::arg("matrix"), py::arg("vector"),
    "Matrix-vector multiply: out = matrix * vector. 'format' is 'R' (row-major) or 'C' (column-major). Returns the new output vector.");

    // void diagmat_vec(const f64* D, bool invD, const f64* x,
    //                  int size, f64* out) -> out
    linalg_m.def("diagmat_vec", [](
        py::array_t<const f64, py::array::c_style | py::array::forcecast> diag_arr,
        bool invD,
        py::array_t<const f64, py::array::c_style | py::array::forcecast> x_arr)
    {
        if (diag_arr.size() != x_arr.size()) throw std::runtime_error("Matrix array size must be same as vector size.");

        py::array_t<f64> out_arr(x_arr.size());
        Linalg::diagmat_vec(
            diag_arr.data(),
            invD,
            x_arr.data(),
            x_arr.size(),
            out_arr.mutable_data()
        );
        return out_arr;
    },
    py::arg("diag"), py::arg("invD"), py::arg("x"),
    "Diagonal matrix-vector multiply: out = diag * x. 'invD' indicates whether to use the inverse of the diagonal.");

    // void dsaxpy(int size, const f64* x, const f64* y,
    //             const f64* D, f64 beta, bool invD, f64* out) -> out
    linalg_m.def("dsaxpy", [](
        py::array_t<const f64, py::array::c_style | py::array::forcecast> x_arr,
        py::array_t<const f64, py::array::c_style | py::array::forcecast> y_arr,
        py::array_t<const f64, py::array::c_style | py::array::forcecast> D_arr,
        f64 beta,
        bool invD)
    {
        if (D_arr.size() != x_arr.size() && D_arr.size() != y_arr.size()) throw std::runtime_error("Matrix array size must be same as vectors sizes.");

        py::array_t<f64> out_arr(x_arr.size());
        Linalg::dsaxpy(
            x_arr.size(),
            x_arr.data(),
            y_arr.data(),
            D_arr.data(),
            beta,
            invD,
            out_arr.mutable_data()
        );
        return out_arr;
    },
    py::arg("x"), py::arg("y"), py::arg("D"), py::arg("beta"), py::arg("invD"),
    "Perform diagonal scaled axpy: out = D * (x + beta * y) or out = D^(-1) * (x + beta * y).");

    // void dgmv(int size, const f64* x, const f64* y,
    //           const f64* D, f64 beta, bool invD, f64* out) -> out
    linalg_m.def("dgmv", [](
        py::array_t<const f64, py::array::c_style | py::array::forcecast> x_arr,
        py::array_t<const f64, py::array::c_style | py::array::forcecast> y_arr,
        py::array_t<const f64, py::array::c_style | py::array::forcecast> D_arr,
        f64 beta,
        bool invD)
    {
        if (D_arr.size() != x_arr.size() && D_arr.size() != y_arr.size()) throw std::runtime_error("Matrix array size must be same as vectors sizes.");

        py::array_t<f64> out_arr(x_arr.size());
        Linalg::dgmv(
            x_arr.size(),
            x_arr.data(),
            y_arr.data(),
            D_arr.data(),
            beta,
            invD,
            out_arr.mutable_data()
        );
        return out_arr;
    },
    py::arg("x"), py::arg("y"), py::arg("D"), py::arg("beta"), py::arg("invD"),
    "Perform diagonal general matrix vector: out = D * x + beta * y or D^(-1) * x + beta * y.");
 
}
