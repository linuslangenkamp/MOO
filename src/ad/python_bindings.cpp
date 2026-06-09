// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2026 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <ad.h>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

struct PyFunction {
    ad::GraphFunction fn;
};

struct PyGraphScalar {
    std::shared_ptr<ad::GraphFunctionBuilder> builder;
    ad::Expr expr;

    PyGraphScalar() = default;
    PyGraphScalar(std::shared_ptr<ad::GraphFunctionBuilder> builder_in, ad::Expr expr_in)
        : builder(std::move(builder_in)),
          expr(expr_in) {}
};

struct PyGraphVector {
    std::shared_ptr<ad::GraphFunctionBuilder> builder;
    ad::VectorExpr expr;

    PyGraphVector() = default;
    PyGraphVector(std::shared_ptr<ad::GraphFunctionBuilder> builder_in, ad::VectorExpr expr_in)
        : builder(std::move(builder_in)),
          expr(expr_in) {}
};

using NamedValues = std::unordered_map<std::string, std::vector<double>>;
using SparsityPairs = std::vector<std::pair<int, int>>;

struct PyPreparedFunction {
    std::shared_ptr<ad::StagedVM> vm;
    ad::StagedVM::Prepared prepared;
};

class PyGraphFunctionBuilder {
public:
    PyGraphFunctionBuilder()
        : builder(std::make_shared<ad::GraphFunctionBuilder>()) {}

    PyGraphScalar constant(double value) {
        ensure_alive();
        return PyGraphScalar(builder, builder->constant(value));
    }

    PyGraphVector inputs(const std::string &name, int size) {
        ensure_alive();
        return PyGraphVector(builder, builder->inputs(name, size));
    }

    PyGraphVector params(const std::string &name, int size) {
        ensure_alive();
        return PyGraphVector(builder, builder->params(name, size));
    }

    PyGraphVector vector(const std::vector<PyGraphScalar> &values) {
        ensure_alive();
        return PyGraphVector(builder, builder->vector(unwrap_scalars(values)));
    }

    PyGraphVector dense_matvec(const std::vector<std::vector<double>> &matrix, const PyGraphVector &rhs) {
        ensure_alive();
        require_same(rhs);
        return PyGraphVector(builder, builder->dense_matvec(make_dense_matrix(matrix), rhs.expr));
    }

    PyGraphVector dense_matvec_values(const std::vector<std::vector<double>> &matrix, const std::vector<PyGraphScalar> &rhs) {
        ensure_alive();
        auto rhs_vector = builder->vector(unwrap_scalars(rhs));
        return PyGraphVector(builder, builder->dense_matvec(make_dense_matrix(matrix), rhs_vector));
    }

    PyGraphVector sparse_matvec(const std::vector<int> &rows, const std::vector<int> &cols, const std::vector<double> &values, const std::pair<int, int> &shape, const PyGraphVector &rhs) {
        ensure_alive();
        require_same(rhs);
        return PyGraphVector(builder, builder->sparse_matvec(ad::SparseMatrix(shape.first, shape.second, rows, cols, values), rhs.expr));
    }

    PyGraphVector sparse_matvec_values(const std::vector<int> &rows,
                                       const std::vector<int> &cols,
                                       const std::vector<double> &values,
                                       const std::pair<int, int> &shape,
                                       const std::vector<PyGraphScalar> &rhs) {
        ensure_alive();
        auto rhs_vector = builder->vector(unwrap_scalars(rhs));
        return PyGraphVector(builder, builder->sparse_matvec(ad::SparseMatrix(shape.first, shape.second, rows, cols, values), rhs_vector));
    }

    PyGraphVector kron_eye_matvec(const std::vector<std::vector<double>> &base, int eye_size, const PyGraphVector &rhs) {
        ensure_alive();
        require_same(rhs);
        return PyGraphVector(builder, builder->kron_eye_matvec(make_dense_matrix(base), eye_size, rhs.expr));
    }

    PyGraphVector kron_eye_matvec_values(const std::vector<std::vector<double>> &base, int eye_size, const std::vector<PyGraphScalar> &rhs) {
        ensure_alive();
        auto rhs_vector = builder->vector(unwrap_scalars(rhs));
        return PyGraphVector(builder, builder->kron_eye_matvec(make_dense_matrix(base), eye_size, rhs_vector));
    }

    PyGraphVector vector_add(const PyGraphVector &lhs, const PyGraphVector &rhs) {
        ensure_alive();
        require_same(lhs);
        require_same(rhs);
        return PyGraphVector(builder, builder->add(lhs.expr, rhs.expr));
    }

    PyGraphVector vector_sub(const PyGraphVector &lhs, const PyGraphVector &rhs) {
        ensure_alive();
        require_same(lhs);
        require_same(rhs);
        return PyGraphVector(builder, builder->sub(lhs.expr, rhs.expr));
    }

    PyGraphVector vector_scale(double factor, const PyGraphVector &rhs) {
        ensure_alive();
        require_same(rhs);
        return PyGraphVector(builder, builder->scale(factor, rhs.expr));
    }

    PyGraphVector concat(const PyGraphVector &lhs, const PyGraphVector &rhs) {
        ensure_alive();
        require_same(lhs);
        require_same(rhs);
        return PyGraphVector(builder, builder->concat(lhs.expr, rhs.expr));
    }

    PyGraphVector slice(const PyGraphVector &rhs, int start, int length, int stride = 1) {
        ensure_alive();
        require_same(rhs);
        return PyGraphVector(builder, builder->slice(rhs.expr, start, length, stride));
    }

    PyGraphScalar at(const PyGraphVector &rhs, int index) {
        ensure_alive();
        require_same(rhs);
        return PyGraphScalar(builder, builder->at(rhs.expr, index));
    }

    PyFunction function(const PyGraphVector &inputs, const PyGraphVector &outputs, const std::vector<PyGraphVector> &params = {}) {
        ensure_alive();
        require_same(inputs);
        require_same(outputs);
        std::vector<ad::VectorExpr> native_params;
        native_params.reserve(params.size());
        for (const auto &param : params) {
            require_same(param);
            native_params.push_back(param.expr);
        }
        PyFunction out;
        out.fn = builder->function(inputs.expr, outputs.expr, native_params);
        builder.reset();
        return out;
    }

private:
    std::shared_ptr<ad::GraphFunctionBuilder> builder;

    void ensure_alive() const {
        if (!builder) {
            throw std::runtime_error("GraphFunctionBuilder has already been frozen into a Function");
        }
    }

    void require_same(const PyGraphScalar &expr) const {
        ensure_alive();
        if (expr.builder.get() != builder.get()) {
            throw std::runtime_error("scalar expression belongs to a different graph function builder");
        }
    }

    void require_same(const PyGraphVector &expr) const {
        ensure_alive();
        if (expr.builder.get() != builder.get()) {
            throw std::runtime_error("vector expression belongs to a different graph function builder");
        }
    }

    std::vector<ad::Expr> unwrap_scalars(const std::vector<PyGraphScalar> &values) const {
        std::vector<ad::Expr> out;
        out.reserve(values.size());
        for (const auto &value : values) {
            require_same(value);
            out.push_back(value.expr);
        }
        return out;
    }

    static ad::DenseMatrix make_dense_matrix(const std::vector<std::vector<double>> &matrix) {
        int rows = static_cast<int>(matrix.size());
        int cols = rows == 0 ? 0 : static_cast<int>(matrix.front().size());
        std::vector<double> flat;
        flat.reserve(static_cast<std::size_t>(rows * std::max(cols, 0)));
        for (const auto &row : matrix) {
            if (static_cast<int>(row.size()) != cols) {
                throw std::runtime_error("dense matrix rows must have equal length");
            }
            flat.insert(flat.end(), row.begin(), row.end());
        }
        return ad::DenseMatrix(rows, cols, std::move(flat));
    }
};

PyGraphScalar require_graph_scalar(py::handle value, const PyGraphScalar &reference) {
    if (py::isinstance<PyGraphScalar>(value)) {
        return value.cast<PyGraphScalar>();
    }
    if (py::isinstance<py::float_>(value) || py::isinstance<py::int_>(value)) {
        return PyGraphScalar(reference.builder, reference.builder->constant(py::cast<double>(value)));
    }
    throw py::type_error("expected AD graph scalar or numeric constant");
}

PyGraphScalar graph_binary(py::handle lhs, py::handle rhs, ad::Op op) {
    PyGraphScalar a;
    PyGraphScalar b;
    if (py::isinstance<PyGraphScalar>(lhs)) {
        a = lhs.cast<PyGraphScalar>();
        b = require_graph_scalar(rhs, a);
    } else if (py::isinstance<PyGraphScalar>(rhs)) {
        b = rhs.cast<PyGraphScalar>();
        a = require_graph_scalar(lhs, b);
    } else {
        throw py::type_error("at least one operand must be an AD graph scalar");
    }
    if (a.builder.get() != b.builder.get()) {
        throw std::runtime_error("graph scalar expressions belong to different graph function builders");
    }
    return PyGraphScalar(a.builder, a.expr.g->binary(op, a.expr, b.expr));
}

PyGraphScalar graph_unary(const PyGraphScalar &x, ad::Op op) {
    return PyGraphScalar(x.builder, x.expr.g->unary(op, x.expr));
}

std::vector<std::pair<int, int>> structural_pairs(const PyFunction &fn, const std::string &wrt) {
    return ad::structural_sparsity(fn.fn, wrt).to_pairs();
}

void add_env_values(ad::EvalEnv &env, const PyFunction &fn, const NamedValues &inputs, const NamedValues &params) {
    for (const auto &[name, values] : inputs) {
        int expected = fn.fn.input_group_size(name);
        if (expected < 0) {
            throw std::runtime_error("unknown input group: " + name);
        }
        if (static_cast<int>(values.size()) != expected) {
            throw std::runtime_error("input group '" + name + "' expects " + std::to_string(expected) + " values, got " + std::to_string(values.size()));
        }
        env.input(name, values.data());
    }
    for (const auto &[name, values] : params) {
        int expected = fn.fn.param_group_size(name);
        if (expected < 0) {
            throw std::runtime_error("unknown parameter group: " + name);
        }
        if (static_cast<int>(values.size()) != expected) {
            throw std::runtime_error("parameter group '" + name + "' expects " + std::to_string(expected) + " values, got " + std::to_string(values.size()));
        }
        env.param(name, values.data());
    }
}

void require_all_groups_bound(const PyFunction &fn, const NamedValues &inputs, const NamedValues &params, const std::string &staged_direction = "") {
    for (const auto &[name, size] : fn.fn.input_groups) {
        (void)size;
        if (inputs.find(name) == inputs.end()) {
            throw std::runtime_error("missing input group: " + name);
        }
    }
    for (const auto &[name, size] : fn.fn.param_groups) {
        (void)size;
        if (name == staged_direction) {
            continue;
        }
        if (params.find(name) == params.end()) {
            throw std::runtime_error("missing parameter group: " + name);
        }
    }
}

std::vector<double> evaluate_vm(const PyFunction &fn, const NamedValues &inputs, const NamedValues &params) {
    require_all_groups_bound(fn, inputs, params);
    ad::EvalEnv env;
    add_env_values(env, fn, inputs, params);
    ad::VM vm(fn.fn);
    std::vector<double> out(fn.fn.output_size(), 0.0);
    vm.evaluate(env, out.data());
    return out;
}

} // namespace

PYBIND11_MODULE(_ad, m) {
    m.doc() = "MOO native AD bindings";

    py::class_<PyGraphScalar>(m, "GraphScalar")
        .def(
            "__add__",
            [](py::handle a, py::handle b) { return graph_binary(a, b, ad::Op::Add); },
            py::is_operator())
        .def(
            "__radd__",
            [](py::handle a, py::handle b) { return graph_binary(b, a, ad::Op::Add); },
            py::is_operator())
        .def(
            "__sub__",
            [](py::handle a, py::handle b) { return graph_binary(a, b, ad::Op::Sub); },
            py::is_operator())
        .def(
            "__rsub__",
            [](py::handle a, py::handle b) { return graph_binary(b, a, ad::Op::Sub); },
            py::is_operator())
        .def(
            "__mul__",
            [](py::handle a, py::handle b) { return graph_binary(a, b, ad::Op::Mul); },
            py::is_operator())
        .def(
            "__rmul__",
            [](py::handle a, py::handle b) { return graph_binary(b, a, ad::Op::Mul); },
            py::is_operator())
        .def(
            "__truediv__",
            [](py::handle a, py::handle b) { return graph_binary(a, b, ad::Op::Div); },
            py::is_operator())
        .def(
            "__rtruediv__",
            [](py::handle a, py::handle b) { return graph_binary(b, a, ad::Op::Div); },
            py::is_operator())
        .def(
            "__neg__",
            [](const PyGraphScalar &x) { return graph_unary(x, ad::Op::Neg); },
            py::is_operator())
        .def("pow_const", [](const PyGraphScalar &x, double p) { return PyGraphScalar(x.builder, x.expr.g->pow_const(x.expr, p)); });

    py::class_<PyGraphVector>(m, "GraphVector")
        .def("__len__", [](const PyGraphVector &v) { return v.expr.size; })
        .def("__getitem__",
             [](const PyGraphVector &v, int i) {
                 if (!v.builder) {
                     throw std::runtime_error("invalid graph vector");
                 }
                 if (i < 0) {
                     i += v.expr.size;
                 }
                 return PyGraphScalar(v.builder, v.builder->at(v.expr, i));
             })
        .def_property_readonly("values",
                               [](const PyGraphVector &v) {
                                   if (!v.builder) {
                                       throw std::runtime_error("invalid graph vector");
                                   }
                                   std::vector<PyGraphScalar> out;
                                   out.reserve(static_cast<std::size_t>(v.expr.size));
                                   for (int i = 0; i < v.expr.size; ++i) {
                                       out.emplace_back(v.builder, v.builder->at(v.expr, i));
                                   }
                                   return out;
                               })
        .def(
            "__add__",
            [](const PyGraphVector &lhs, const PyGraphVector &rhs) {
                if (lhs.builder.get() != rhs.builder.get()) {
                    throw std::runtime_error("graph vectors belong to different graph function builders");
                }
                return PyGraphVector(lhs.builder, lhs.builder->add(lhs.expr, rhs.expr));
            },
            py::is_operator())
        .def(
            "__sub__",
            [](const PyGraphVector &lhs, const PyGraphVector &rhs) {
                if (lhs.builder.get() != rhs.builder.get()) {
                    throw std::runtime_error("graph vectors belong to different graph function builders");
                }
                return PyGraphVector(lhs.builder, lhs.builder->sub(lhs.expr, rhs.expr));
            },
            py::is_operator())
        .def(
            "__mul__",
            [](const PyGraphVector &rhs, double factor) { return PyGraphVector(rhs.builder, rhs.builder->scale(factor, rhs.expr)); },
            py::is_operator())
        .def(
            "__rmul__",
            [](const PyGraphVector &rhs, double factor) { return PyGraphVector(rhs.builder, rhs.builder->scale(factor, rhs.expr)); },
            py::is_operator());

    py::class_<PyGraphFunctionBuilder>(m, "GraphFunctionBuilder")
        .def(py::init<>())
        .def("constant", &PyGraphFunctionBuilder::constant)
        .def("inputs", &PyGraphFunctionBuilder::inputs)
        .def("params", &PyGraphFunctionBuilder::params)
        .def("vector", &PyGraphFunctionBuilder::vector)
        .def("dense_matvec", &PyGraphFunctionBuilder::dense_matvec)
        .def("dense_matvec_values", &PyGraphFunctionBuilder::dense_matvec_values)
        .def("sparse_matvec", &PyGraphFunctionBuilder::sparse_matvec)
        .def("sparse_matvec_values", &PyGraphFunctionBuilder::sparse_matvec_values)
        .def("kron_eye_matvec", &PyGraphFunctionBuilder::kron_eye_matvec)
        .def("kron_eye_matvec_values", &PyGraphFunctionBuilder::kron_eye_matvec_values)
        .def("vector_add", &PyGraphFunctionBuilder::vector_add)
        .def("vector_sub", &PyGraphFunctionBuilder::vector_sub)
        .def("vector_scale", &PyGraphFunctionBuilder::vector_scale)
        .def("concat", &PyGraphFunctionBuilder::concat)
        .def("slice", &PyGraphFunctionBuilder::slice, py::arg("rhs"), py::arg("start"), py::arg("length"), py::arg("stride") = 1)
        .def("at", &PyGraphFunctionBuilder::at)
        .def("function", &PyGraphFunctionBuilder::function, py::arg("inputs"), py::arg("outputs"), py::arg("params") = std::vector<PyGraphVector>{});

    py::class_<PyFunction>(m, "Function")
        .def_property_readonly("input_size", [](const PyFunction &f) { return f.fn.input_size(); })
        .def_property_readonly("param_size", [](const PyFunction &f) { return f.fn.param_size(); })
        .def_property_readonly("output_size", [](const PyFunction &f) { return f.fn.output_size(); })
        .def_property_readonly("has_vector_structure", [](const PyFunction &f) { return f.fn.has_vector_structure(); })
        .def_property_readonly("vector_node_count", [](const PyFunction &f) { return f.fn.vector_node_count(); })
        .def("evaluate", &evaluate_vm, py::arg("inputs") = NamedValues{}, py::arg("params") = NamedValues{})
        .def(
            "prepare",
            [](const PyFunction &f, const NamedValues &inputs, const NamedValues &params, const std::string &direction) {
                require_all_groups_bound(f, inputs, params, direction);
                ad::EvalEnv env;
                add_env_values(env, f, inputs, params);
                PyPreparedFunction prepared;
                prepared.vm = std::make_shared<ad::StagedVM>(f.fn, direction);
                prepared.prepared = prepared.vm->prepare(env);
                return prepared;
            },
            py::arg("inputs") = NamedValues{},
            py::arg("params") = NamedValues{},
            py::arg("direction") = "v")
        .def(
            "staged_apply",
            [](const PyFunction &f, const std::vector<double> &direction_values, const NamedValues &inputs, const NamedValues &params, const std::string &direction) {
                PyPreparedFunction prepared;
                require_all_groups_bound(f, inputs, params, direction);
                ad::EvalEnv env;
                add_env_values(env, f, inputs, params);
                prepared.vm = std::make_shared<ad::StagedVM>(f.fn, direction);
                int expected = f.fn.param_group_size(direction);
                if (expected < 0) {
                    throw std::runtime_error("unknown direction parameter group: " + direction);
                }
                if (static_cast<int>(direction_values.size()) != expected) {
                    throw std::runtime_error("direction group '" + direction + "' expects " + std::to_string(expected) + " values, got " + std::to_string(direction_values.size()));
                }
                prepared.prepared = prepared.vm->prepare(env);
                std::vector<double> out(f.fn.output_size(), 0.0);
                prepared.prepared.apply(direction_values.data(), out.data());
                return out;
            },
            py::arg("direction_values"),
            py::arg("inputs") = NamedValues{},
            py::arg("params") = NamedValues{},
            py::arg("direction") = "v")
        .def(
            "forward_diff",
            [](const PyFunction &f, const std::string &wrt, const std::string &direction) { return PyFunction{ad::forward_diff(f.fn, wrt, direction)}; },
            py::arg("wrt") = "x",
            py::arg("direction") = "v")
        .def(
            "reverse_diff",
            [](const PyFunction &f, const std::string &lambda, const std::string &wrt) { return PyFunction{ad::reverse_diff(f.fn, lambda, wrt)}; },
            py::arg("lambda_name") = "lambda",
            py::arg("wrt") = "x")
        .def(
            "exact_derivative_plan",
            [](const PyFunction &f, const std::string &wrt, const std::string &direction, const std::string &lambda_name) {
                auto plan = ad::exact_derivative_plan(f.fn, wrt, direction, lambda_name);
                py::dict out;
                out["jvp"] = PyFunction{std::move(plan.jvp)};
                out["grad"] = PyFunction{std::move(plan.grad)};
                out["hvp"] = PyFunction{std::move(plan.hvp)};
                out["jacobian_sparsity"] = std::move(plan.jacobian_sparsity);
                out["hessian_sparsity"] = std::move(plan.hessian_sparsity);
                out["hessian_full_sparsity"] = std::move(plan.hessian_full_sparsity);
                out["jacobian_colors"] = std::move(plan.jacobian_colors);
                out["hessian_colors"] = std::move(plan.hessian_colors);
                out["jacobian_color_count"] = plan.jacobian_color_count;
                out["hessian_color_count"] = plan.hessian_color_count;
                return out;
            },
            py::arg("wrt") = "x",
            py::arg("direction") = "v",
            py::arg("lambda_name") = "lambda")
        .def("to_c", [](const PyFunction &f, const std::string &name) { return ad::to_c(f.fn, name); })
        .def(
            "to_staged_c",
            [](const PyFunction &f, const std::string &name, const std::string &direction) { return ad::to_staged_c(f.fn, name, direction); },
            py::arg("name"),
            py::arg("direction") = "v")
        .def(
            "to_sparse_coefficients_c",
            [](const PyFunction &f, const std::string &seed_name, const SparsityPairs &pairs, const std::string &name) {
                return ad::to_sparse_coefficients_c(f.fn, seed_name, pairs, name);
            },
            py::arg("seed_name"),
            py::arg("pairs"),
            py::arg("name"))
        .def(
            "to_sparse_jacobian_c",
            [](const PyFunction &f, const std::string &wrt, const SparsityPairs &pairs, const std::string &name) { return ad::to_sparse_jacobian_c(f.fn, wrt, pairs, name); },
            py::arg("wrt"),
            py::arg("pairs"),
            py::arg("name"))
        .def(
            "to_sparse_hessian_c",
            [](const PyFunction &f, const std::string &direction, const SparsityPairs &pairs, const std::string &name) {
                return ad::to_sparse_hessian_c(f.fn, direction, pairs, name);
            },
            py::arg("direction"),
            py::arg("pairs"),
            py::arg("name"))
        .def(
            "emit_exact_derivative_code",
            [](const PyFunction &f,
               const std::string &wrt,
               const std::string &direction,
               const std::string &lambda_name,
               const std::string &value_name,
               const std::string &jvp_name,
               const std::string &hvp_name) {
                auto code = ad::emit_exact_derivative_code(f.fn, wrt, direction, lambda_name, value_name, jvp_name, hvp_name);
                py::dict out;
                out["value"] = std::move(code.value);
                out["jvp"] = std::move(code.jvp);
                out["hvp"] = std::move(code.hvp);
                out["jacobian"] = std::move(code.jacobian);
                out["hessian"] = std::move(code.hessian);
                out["jacobian_sparsity"] = std::move(code.jacobian_sparsity);
                out["hessian_sparsity"] = std::move(code.hessian_sparsity);
                out["jacobian_colors"] = std::move(code.jacobian_colors);
                out["hessian_colors"] = std::move(code.hessian_colors);
                out["jacobian_color_count"] = code.jacobian_color_count;
                out["hessian_color_count"] = code.hessian_color_count;
                out["value_bytes"] = code.value_bytes;
                out["jvp_bytes"] = code.jvp_bytes;
                out["hvp_bytes"] = code.hvp_bytes;
                out["jacobian_bytes"] = code.jacobian_bytes;
                out["hessian_bytes"] = code.hessian_bytes;
                return out;
            },
            py::arg("wrt"),
            py::arg("direction"),
            py::arg("lambda_name"),
            py::arg("value_name"),
            py::arg("jvp_name"),
            py::arg("hvp_name"))
        .def(
            "coloring",
            [](const PyFunction &f, const SparsityPairs &pairs, int column_count) {
                auto colors = ad::greedy_column_coloring(column_count, pairs);
                int color_count = 0;
                for (int color : colors) {
                    color_count = std::max(color_count, color + 1);
                }
                py::dict out;
                out["colors"] = colors;
                out["color_count"] = color_count;
                return out;
            },
            py::arg("pairs"),
            py::arg("column_count"))
        .def("structural_sparsity", &structural_pairs)
        .def(
            "jacobian_sparsity",
            [](const PyFunction &f, const std::string &wrt) { return ad::jacobian_sparsity(f.fn, wrt).to_pairs(); },
            py::arg("wrt") = "x")
        .def(
            "hessian_sparsity",
            [](const PyFunction &f, const std::string &direction) { return ad::hessian_sparsity(f.fn, direction); },
            py::arg("direction") = "v")
        .def(
            "hessian_sparsity_full",
            [](const PyFunction &f, const std::string &direction) { return ad::hessian_sparsity_full(f.fn, direction); },
            py::arg("direction") = "v");

    py::class_<PyPreparedFunction>(m, "PreparedFunction")
        .def(
            "apply",
            [](const PyPreparedFunction &p, const std::vector<double> &direction_values) {
                int expected = ad::param_group_size(p.vm->f, p.vm->direction_name);
                if (static_cast<int>(direction_values.size()) != expected) {
                    throw std::runtime_error("direction group '" + p.vm->direction_name + "' expects " + std::to_string(expected) + " values, got " +
                                             std::to_string(direction_values.size()));
                }
                std::vector<double> out(p.vm->f.output_size(), 0.0);
                p.prepared.apply(direction_values.data(), out.data());
                return out;
            },
            py::arg("direction_values"));

    m.def("derivative_callback_mode", &ad::derivative_callback_mode, py::arg("strategy"), py::arg("pairs"), py::arg("colors"));
    m.def("derivative_section_keys", &ad::derivative_section_keys, py::arg("prefix"), py::arg("jac_mode"), py::arg("hes_mode"));
    m.def("render_jacobian_callback_body",
          &ad::render_jacobian_callback_body,
          py::arg("mode"),
          py::arg("direct_call"),
          py::arg("input_size"),
          py::arg("output_size"),
          py::arg("pairs"),
          py::arg("colors"),
          py::arg("colored_call"));
    m.def("render_hessian_callback_body",
          &ad::render_hessian_callback_body,
          py::arg("mode"),
          py::arg("direct_body"),
          py::arg("input_size"),
          py::arg("tmp_size"),
          py::arg("pairs"),
          py::arg("colors"),
          py::arg("prepare_body"),
          py::arg("apply_call"),
          py::arg("buf_indices") = std::vector<int>{});

    m.def("sin", [](const PyGraphScalar &x) { return graph_unary(x, ad::Op::Sin); });
    m.def("cos", [](const PyGraphScalar &x) { return graph_unary(x, ad::Op::Cos); });
    m.def("tan", [](const PyGraphScalar &x) { return graph_unary(x, ad::Op::Tan); });
    m.def("exp", [](const PyGraphScalar &x) { return graph_unary(x, ad::Op::Exp); });
    m.def("log", [](const PyGraphScalar &x) { return graph_unary(x, ad::Op::Log); });
    m.def("pow_const", [](const PyGraphScalar &x, double p) { return PyGraphScalar(x.builder, x.expr.g->pow_const(x.expr, p)); });
}
