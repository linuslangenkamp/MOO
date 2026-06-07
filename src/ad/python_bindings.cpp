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
// You should have received a copy of the GNU General Public License
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

struct PyExpr {
    std::shared_ptr<ad::Graph> graph;
    ad::Expr expr;

    PyExpr() = default;
    PyExpr(std::shared_ptr<ad::Graph> graph, ad::Expr expr)
        : graph(std::move(graph)), expr(expr)
    {
    }
};

struct PyVector {
    std::string name;
    std::vector<PyExpr> values;
};

struct PyFunction {
    ad::GraphFunction fn;
};

using NamedValues = std::unordered_map<std::string, std::vector<double>>;
using SparsityPairs = std::vector<std::pair<int, int>>;

struct PyPreparedFunction {
    std::shared_ptr<ad::StagedVM> vm;
    ad::StagedVM::Prepared prepared;
};

class PyGraphBuilder {
public:
    PyGraphBuilder()
        : graph(std::make_shared<ad::Graph>())
    {
    }

    PyExpr constant(double value)
    {
        return PyExpr(graph, graph->constant(value));
    }

    PyVector inputs(const std::string& name, int size)
    {
        auto vec = graph->inputs(name, size);
        return wrap_vector(vec);
    }

    PyVector params(const std::string& name, int size)
    {
        auto vec = graph->params(name, size);
        return wrap_vector(vec);
    }

    PyFunction function(const PyVector& inputs, const std::vector<PyExpr>& outputs, const std::vector<PyVector>& params = {})
    {
        ensure_alive();
        ad::NamedVector native_inputs = unwrap_vector(inputs);
        std::vector<ad::Expr> native_outputs;
        native_outputs.reserve(outputs.size());
        for (const auto& out : outputs) {
            require_same(out);
            native_outputs.push_back(out.expr);
        }
        std::vector<ad::NamedVector> native_params;
        native_params.reserve(params.size());
        for (const auto& param : params) {
            native_params.push_back(unwrap_vector(param));
        }

        PyFunction out;
        out.fn = ad::function_from(std::move(*graph), native_inputs, native_outputs, native_params);
        graph.reset();
        return out;
    }

private:
    std::shared_ptr<ad::Graph> graph;

    void ensure_alive() const
    {
        if (!graph) {
            throw std::runtime_error("GraphBuilder has already been frozen into a Function");
        }
    }

    void require_same(const PyExpr& expr) const
    {
        ensure_alive();
        if (expr.graph.get() != graph.get()) {
            throw std::runtime_error("expression belongs to a different graph");
        }
    }

    PyVector wrap_vector(const ad::NamedVector& native)
    {
        PyVector out;
        out.name = native.name;
        out.values.reserve(native.values.size());
        for (const auto& expr : native.values) {
            out.values.emplace_back(graph, expr);
        }
        return out;
    }

    ad::NamedVector unwrap_vector(const PyVector& vec) const
    {
        ad::NamedVector out;
        out.name = vec.name;
        out.values.reserve(vec.values.size());
        for (const auto& expr : vec.values) {
            require_same(expr);
            out.values.push_back(expr.expr);
        }
        return out;
    }
};

PyExpr require_expr(py::handle value, const PyExpr& reference)
{
    if (py::isinstance<PyExpr>(value)) {
        return value.cast<PyExpr>();
    }
    if (py::isinstance<py::float_>(value) || py::isinstance<py::int_>(value)) {
        return PyExpr(reference.graph, reference.graph->constant(py::cast<double>(value)));
    }
    throw py::type_error("expected AD expression or numeric constant");
}

PyExpr binary(py::handle lhs, py::handle rhs, ad::Op op)
{
    PyExpr a;
    PyExpr b;
    if (py::isinstance<PyExpr>(lhs)) {
        a = lhs.cast<PyExpr>();
        b = require_expr(rhs, a);
    } else if (py::isinstance<PyExpr>(rhs)) {
        b = rhs.cast<PyExpr>();
        a = require_expr(lhs, b);
    } else {
        throw py::type_error("at least one operand must be an AD expression");
    }
    if (a.graph.get() != b.graph.get()) {
        throw std::runtime_error("expressions belong to different graphs");
    }
    return PyExpr(a.graph, a.graph->binary(op, a.expr, b.expr));
}

PyExpr unary(const PyExpr& x, ad::Op op)
{
    return PyExpr(x.graph, x.graph->unary(op, x.expr));
}

std::vector<std::pair<int, int>> structural_pairs(const PyFunction& fn, const std::string& wrt)
{
    return ad::structural_sparsity(fn.fn, wrt).to_pairs();
}

int group_size(const std::vector<std::pair<std::string, int>>& groups, const std::string& name)
{
    for (const auto& group : groups) {
        if (group.first == name) {
            return group.second;
        }
    }
    return -1;
}

void add_env_values(ad::EvalEnv& env, const PyFunction& fn, const NamedValues& inputs, const NamedValues& params)
{
    for (const auto& [name, values] : inputs) {
        int expected = group_size(fn.fn.input_groups, name);
        if (expected < 0) {
            throw std::runtime_error("unknown input group: " + name);
        }
        if (static_cast<int>(values.size()) != expected) {
            throw std::runtime_error("input group '" + name + "' expects " + std::to_string(expected)
                                     + " values, got " + std::to_string(values.size()));
        }
        env.input(name, values.data());
    }
    for (const auto& [name, values] : params) {
        int expected = group_size(fn.fn.param_groups, name);
        if (expected < 0) {
            throw std::runtime_error("unknown parameter group: " + name);
        }
        if (static_cast<int>(values.size()) != expected) {
            throw std::runtime_error("parameter group '" + name + "' expects " + std::to_string(expected)
                                     + " values, got " + std::to_string(values.size()));
        }
        env.param(name, values.data());
    }
}

void require_all_groups_bound(const PyFunction& fn, const NamedValues& inputs, const NamedValues& params,
                              const std::string& staged_direction = "")
{
    for (const auto& [name, size] : fn.fn.input_groups) {
        (void)size;
        if (inputs.find(name) == inputs.end()) {
            throw std::runtime_error("missing input group: " + name);
        }
    }
    for (const auto& [name, size] : fn.fn.param_groups) {
        (void)size;
        if (name == staged_direction) {
            continue;
        }
        if (params.find(name) == params.end()) {
            throw std::runtime_error("missing parameter group: " + name);
        }
    }
}

std::vector<double> evaluate_vm(const PyFunction& fn, const NamedValues& inputs, const NamedValues& params)
{
    require_all_groups_bound(fn, inputs, params);
    ad::EvalEnv env;
    add_env_values(env, fn, inputs, params);
    ad::VM vm(fn.fn);
    std::vector<double> out(fn.fn.output_size(), 0.0);
    vm.evaluate(env, out.data());
    return out;
}

} // namespace

PYBIND11_MODULE(_ad, m)
{
    m.doc() = "MOO native AD bindings";

    py::class_<PyExpr>(m, "Expr")
        .def("__add__", [](py::handle a, py::handle b) { return binary(a, b, ad::Op::Add); }, py::is_operator())
        .def("__radd__", [](py::handle a, py::handle b) { return binary(b, a, ad::Op::Add); }, py::is_operator())
        .def("__sub__", [](py::handle a, py::handle b) { return binary(a, b, ad::Op::Sub); }, py::is_operator())
        .def("__rsub__", [](py::handle a, py::handle b) { return binary(b, a, ad::Op::Sub); }, py::is_operator())
        .def("__mul__", [](py::handle a, py::handle b) { return binary(a, b, ad::Op::Mul); }, py::is_operator())
        .def("__rmul__", [](py::handle a, py::handle b) { return binary(b, a, ad::Op::Mul); }, py::is_operator())
        .def("__truediv__", [](py::handle a, py::handle b) { return binary(a, b, ad::Op::Div); }, py::is_operator())
        .def("__rtruediv__", [](py::handle a, py::handle b) { return binary(b, a, ad::Op::Div); }, py::is_operator())
        .def("__neg__", [](const PyExpr& x) { return unary(x, ad::Op::Neg); }, py::is_operator())
        .def("pow_const", [](const PyExpr& x, double p) { return PyExpr(x.graph, x.graph->pow_const(x.expr, p)); });

    py::class_<PyVector>(m, "Vector")
        .def("__len__", [](const PyVector& v) { return v.values.size(); })
        .def("__getitem__", [](const PyVector& v, std::size_t i) {
            if (i >= v.values.size()) {
                throw py::index_error();
            }
            return v.values[i];
        })
        .def_property_readonly("name", [](const PyVector& v) { return v.name; })
        .def_property_readonly("values", [](const PyVector& v) { return v.values; });

    py::class_<PyGraphBuilder>(m, "GraphBuilder")
        .def(py::init<>())
        .def("constant", &PyGraphBuilder::constant)
        .def("inputs", &PyGraphBuilder::inputs)
        .def("params", &PyGraphBuilder::params)
        .def("function", &PyGraphBuilder::function, py::arg("inputs"), py::arg("outputs"), py::arg("params") = std::vector<PyVector>{});

    py::class_<PyFunction>(m, "Function")
        .def_property_readonly("input_size", [](const PyFunction& f) { return f.fn.input_size(); })
        .def_property_readonly("param_size", [](const PyFunction& f) { return f.fn.param_size(); })
        .def_property_readonly("output_size", [](const PyFunction& f) { return f.fn.output_size(); })
        .def("evaluate", &evaluate_vm, py::arg("inputs") = NamedValues{}, py::arg("params") = NamedValues{})
        .def("prepare", [](const PyFunction& f, const NamedValues& inputs, const NamedValues& params, const std::string& direction) {
            require_all_groups_bound(f, inputs, params, direction);
            ad::EvalEnv env;
            add_env_values(env, f, inputs, params);
            PyPreparedFunction prepared;
            prepared.vm = std::make_shared<ad::StagedVM>(f.fn, direction);
            prepared.prepared = prepared.vm->prepare(env);
            return prepared;
        }, py::arg("inputs") = NamedValues{}, py::arg("params") = NamedValues{}, py::arg("direction") = "v")
        .def("staged_apply", [](const PyFunction& f, const std::vector<double>& direction_values,
                                const NamedValues& inputs, const NamedValues& params, const std::string& direction) {
            PyPreparedFunction prepared;
            require_all_groups_bound(f, inputs, params, direction);
            ad::EvalEnv env;
            add_env_values(env, f, inputs, params);
            prepared.vm = std::make_shared<ad::StagedVM>(f.fn, direction);
            int expected = group_size(f.fn.param_groups, direction);
            if (expected < 0) {
                throw std::runtime_error("unknown direction parameter group: " + direction);
            }
            if (static_cast<int>(direction_values.size()) != expected) {
                throw std::runtime_error("direction group '" + direction + "' expects " + std::to_string(expected)
                                         + " values, got " + std::to_string(direction_values.size()));
            }
            prepared.prepared = prepared.vm->prepare(env);
            std::vector<double> out(f.fn.output_size(), 0.0);
            prepared.prepared.apply(direction_values.data(), out.data());
            return out;
        }, py::arg("direction_values"), py::arg("inputs") = NamedValues{}, py::arg("params") = NamedValues{}, py::arg("direction") = "v")
        .def("forward_diff", [](const PyFunction& f, const std::string& wrt, const std::string& direction) {
            return PyFunction{ad::forward_diff(f.fn, wrt, direction)};
        }, py::arg("wrt") = "x", py::arg("direction") = "v")
        .def("reverse_diff", [](const PyFunction& f, const std::string& lambda, const std::string& wrt) {
            return PyFunction{ad::reverse_diff(f.fn, lambda, wrt)};
        }, py::arg("lambda_name") = "lambda", py::arg("wrt") = "x")
        .def("to_c", [](const PyFunction& f, const std::string& name) { return ad::to_c(f.fn, name); })
        .def("to_staged_c", [](const PyFunction& f, const std::string& name, const std::string& direction) {
            return ad::to_staged_c(f.fn, name, direction);
        }, py::arg("name"), py::arg("direction") = "v")
        .def("to_sparse_coefficients_c", [](const PyFunction& f, const std::string& seed_name, const SparsityPairs& pairs, const std::string& name) {
            return ad::to_sparse_coefficients_c(f.fn, seed_name, pairs, name);
        }, py::arg("seed_name"), py::arg("pairs"), py::arg("name"))
        .def("to_sparse_jacobian_c", [](const PyFunction& f, const std::string& wrt, const SparsityPairs& pairs, const std::string& name) {
            return ad::to_sparse_jacobian_c(f.fn, wrt, pairs, name);
        }, py::arg("wrt"), py::arg("pairs"), py::arg("name"))
        .def("to_sparse_hessian_c", [](const PyFunction& f, const std::string& direction, const SparsityPairs& pairs, const std::string& name) {
            return ad::to_sparse_hessian_c(f.fn, direction, pairs, name);
        }, py::arg("direction"), py::arg("pairs"), py::arg("name"))
        .def("coloring", [](const PyFunction& f, const SparsityPairs& pairs, int column_count) {
            auto colors = ad::greedy_column_coloring(column_count, pairs);
            int color_count = 0;
            for (int color : colors) {
                color_count = std::max(color_count, color + 1);
            }
            py::dict out;
            out["colors"] = colors;
            out["color_count"] = color_count;
            return out;
        }, py::arg("pairs"), py::arg("column_count"))
        .def("structural_sparsity", &structural_pairs)
        .def("jacobian_sparsity", [](const PyFunction& f, const std::string& wrt) {
            return ad::jacobian_sparsity(f.fn, wrt).to_pairs();
        }, py::arg("wrt") = "x")
        .def("hessian_sparsity", [](const PyFunction& f, const std::string& direction) {
            return ad::hessian_sparsity(f.fn, direction);
        }, py::arg("direction") = "v")
        .def("hessian_sparsity_full", [](const PyFunction& f, const std::string& direction) {
            return ad::hessian_sparsity_full(f.fn, direction);
        }, py::arg("direction") = "v");

    py::class_<PyPreparedFunction>(m, "PreparedFunction")
        .def("apply", [](const PyPreparedFunction& p, const std::vector<double>& direction_values) {
            int expected = ad::param_group_size(p.vm->f, p.vm->direction_name);
            if (static_cast<int>(direction_values.size()) != expected) {
                throw std::runtime_error("direction group '" + p.vm->direction_name + "' expects " + std::to_string(expected)
                                         + " values, got " + std::to_string(direction_values.size()));
            }
            std::vector<double> out(p.vm->f.output_size(), 0.0);
            p.prepared.apply(direction_values.data(), out.data());
            return out;
        }, py::arg("direction_values"));

    m.def("sin", [](const PyExpr& x) { return unary(x, ad::Op::Sin); });
    m.def("cos", [](const PyExpr& x) { return unary(x, ad::Op::Cos); });
    m.def("tan", [](const PyExpr& x) { return unary(x, ad::Op::Tan); });
    m.def("exp", [](const PyExpr& x) { return unary(x, ad::Op::Exp); });
    m.def("log", [](const PyExpr& x) { return unary(x, ad::Op::Log); });
    m.def("pow_const", [](const PyExpr& x, double p) { return PyExpr(x.graph, x.graph->pow_const(x.expr, p)); });
}
