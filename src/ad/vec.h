// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_VEC_H
#define MOO_AD_VEC_H

#include "expr.h"
#include "matrix.h"

#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

namespace ad {

namespace detail {
struct VecNode;
const std::shared_ptr<const VecNode> &vec_node(const Vec &vec);
Vec vec_from_node(std::shared_ptr<const VecNode> node);
} // namespace detail

class Vec {
public:
    Vec() = default;
    explicit Vec(std::vector<Expr> elements);
    Vec(std::initializer_list<Expr> elements);

    bool valid() const;
    bool empty() const;
    int size() const;
    NodeId node_id() const;
    GraphNodeKind node_kind() const;

    Expr operator[](int index) const;
    Vec slice(int start, int length) const;
    Vars variables() const;
    Params parameters() const;
    void eval(const Values &values, EvalWorkspace &workspace, double *out, int out_size) const;
    Vec forward_diff(const Vars &wrt, const Vec &seed) const;
    Vec reverse_diff(const Vars &wrt, const Vec &lambda) const;

private:
    explicit Vec(std::shared_ptr<const detail::VecNode> node);

    std::shared_ptr<const detail::VecNode> node_;

    friend const std::shared_ptr<const detail::VecNode> &detail::vec_node(const Vec &vec);
    friend Vec detail::vec_from_node(std::shared_ptr<const detail::VecNode> node);
};

Vec vec(std::vector<Expr> elements);
Vec vec_variable(const std::string &label, int size);
Vec vec_parameter(const std::string &label, int size);
Vec vec_constant(const std::vector<double> &values);

Vec concat(const Vec &lhs, const Vec &rhs);
Vec gather(const Vec &source, std::vector<int> indices);
Vec scatter_add(const Vec &values, std::vector<int> indices, int output_size);
Expr sum(const Vec &values);
Expr dot(const Vec &lhs, const Vec &rhs);

Vec operator+(const Vec &lhs, const Vec &rhs);
Vec operator-(const Vec &lhs, const Vec &rhs);
Vec operator*(const Vec &lhs, const Vec &rhs);
Vec operator/(const Vec &lhs, const Vec &rhs);

Vec operator+(const Vec &lhs, double rhs);
Vec operator-(const Vec &lhs, double rhs);
Vec operator*(const Vec &lhs, double rhs);
Vec operator/(const Vec &lhs, double rhs);
Vec operator+(double lhs, const Vec &rhs);
Vec operator-(double lhs, const Vec &rhs);
Vec operator*(double lhs, const Vec &rhs);
Vec operator/(double lhs, const Vec &rhs);

Vec operator*(const Expr &lhs, const Vec &rhs);
Vec operator*(const Vec &lhs, const Expr &rhs);

Vec operator*(const DenseMatrix &matrix, const Vec &rhs);
Vec operator*(const SparseMatrix &matrix, const Vec &rhs);

Vec sin(const Vec &expr);
Vec cos(const Vec &expr);
Vec tan(const Vec &expr);
Vec exp(const Vec &expr);
Vec log(const Vec &expr);
Vec sigmoid(const Vec &expr);

} // namespace ad

#endif // MOO_AD_VEC_H
