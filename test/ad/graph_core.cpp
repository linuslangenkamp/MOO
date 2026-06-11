// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool check(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        return false;
    }
    return true;
}

bool throws_cleanly(void (*fn)()) {
    try {
        fn();
    } catch (const std::exception &) {
        return true;
    }
    return false;
}

bool has_child_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return child.kind == kind;
    });
}

std::string read_file(const std::string &path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open " + path);
    }
    return std::string(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

std::string lower_ascii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool test_scalar_symbol_identity() {
    ad::Expr x1 = ad::variable("x");
    ad::Expr x2 = ad::variable("x");
    ad::Expr p = ad::parameter("x");

    bool ok = true;
    ok &= check(x1.var().label() == x2.var().label(), "same-name scalar variables should keep labels");
    ok &= check(x1.var().id() != x2.var().id(), "same-name scalar variables must have different IDs");
    ok &= check(x1.is_variable() && !x1.is_parameter(), "scalar variable kind is wrong");
    ok &= check(p.is_parameter() && !p.is_variable(), "scalar parameter kind is wrong");
    ok &= check(p.param().label() == x1.var().label(), "parameter label should be independent diagnostic text");
    ok &= check(ad::inspect(x1).kind == ad::GraphNodeKind::ScalarVariable, "scalar variable inspection kind is wrong");
    ok &= check(ad::inspect(p).kind == ad::GraphNodeKind::ScalarParameter, "scalar parameter inspection kind is wrong");
    return ok;
}

bool test_vector_symbol_identity() {
    ad::Vec X1 = ad::vec_variable("X", 3);
    ad::Vec X2 = ad::vec_variable("X", 3);
    ad::Vars v1 = X1.variables();
    ad::Vars v2 = X2.variables();

    bool ok = true;
    ok &= check(v1.size() == 3 && v2.size() == 3, "vector variables should expose component IDs");
    ok &= check(v1[0].label() == "X" && v2[0].label() == "X", "vector variable labels should be preserved");
    ok &= check(v1[0].component() == 0 && v1[2].component() == 2, "vector variable component indices are wrong");
    ok &= check(v1[0].id() != v2[0].id(), "same-name vector variables must have different component IDs");
    ok &= check(v1[0].group_id() != v2[0].group_id(), "same-name vector variables must have different group IDs");

    ad::Vec P = ad::vec_parameter("P", 4);
    ad::Params params = P.parameters();
    ok &= check(params.size() == 4, "vector parameter should expose component IDs");
    ok &= check(params[0].label() == "P" && params[3].component() == 3, "vector parameter metadata is wrong");
    ok &= check(ad::inspect(P).kind == ad::GraphNodeKind::VectorParameter, "vector parameter inspection kind is wrong");
    return ok;
}

bool test_matrix_validation() {
    bool ok = true;
    ok &= check(throws_cleanly([] { (void)ad::DenseMatrix(2, 3, {1.0, 2.0}); }), "bad dense matrix dimensions should throw");
    ok &= check(throws_cleanly([] { (void)ad::DenseMatrix(-1, 0, {}); }), "negative dense matrix dimensions should throw");
    ok &= check(throws_cleanly([] { (void)ad::SparseMatrix(2, 2, {0, 2}, {0, 1}, {1.0, 2.0}); }), "out-of-bounds sparse entry should throw");
    ok &= check(throws_cleanly([] { (void)ad::SparseMatrix(2, 2, {0}, {0, 1}, {1.0}); }), "mismatched sparse triplet arrays should throw");
    ok &= check(ad::DenseMatrix(2, 2, {1.0, 0.0, 0.0, 1.0})(1, 1) == 1.0, "dense matrix indexing is wrong");
    ok &= check(ad::SparseMatrix(2, 3, {0, 1}, {2, 0}, {5.0, -1.0}).nnz() == 2, "sparse matrix nnz is wrong");
    return ok;
}

bool test_structured_matvec_and_unary_graph() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec B = ad::vec_variable("B", 2);
    ad::DenseMatrix A(2, 3, {1.0, 0.0, -2.0,
                             0.0, 3.0, 0.5});

    ad::Vec AX = A * X;
    ad::GraphInfo ax = ad::inspect(AX);

    bool ok = true;
    ok &= check(ax.kind == ad::GraphNodeKind::DenseMatVec, "D * X should create a DenseMatVec node");
    ok &= check(ax.size == 2, "DenseMatVec output size is wrong");
    ok &= check(has_child_kind(ax, ad::GraphNodeKind::VectorVariable), "DenseMatVec should keep the vector variable child");
    ok &= check(ad::count_nodes(ax, ad::GraphNodeKind::VectorElement) == 0, "DenseMatVec should not scalar-lower to vector elements");

    ad::Vec G = ad::sigmoid(AX + B);
    ad::GraphInfo g = ad::inspect(G);
    ok &= check(g.kind == ad::GraphNodeKind::VectorUnary && g.op == "sigmoid", "sigmoid should create a vector unary node");
    ok &= check(ad::count_nodes(g, ad::GraphNodeKind::DenseMatVec) == 1, "sigmoid(A * X + B) should retain DenseMatVec");
    ok &= check(ad::count_nodes(g, ad::GraphNodeKind::VectorBinary) == 1, "sigmoid(A * X + B) should retain vector add");
    ok &= check(ad::count_nodes(g, ad::GraphNodeKind::VectorElement) == 0, "sigmoid(A * X + B) should not scalar-lower");

    ad::SparseMatrix S(2, 3, {0, 1}, {0, 2}, {2.0, -1.0});
    ad::GraphInfo sx = ad::inspect(S * X);
    ok &= check(sx.kind == ad::GraphNodeKind::SparseMatVec, "S * X should create a SparseMatVec node");
    return ok;
}

bool test_slice_concat_scale_and_reductions() {
    ad::Vec X = ad::vec_variable("X", 5);
    ad::Expr h = ad::parameter("h");

    ad::Vec scaled = h * X;
    ad::GraphInfo scaled_info = ad::inspect(scaled);

    bool ok = true;
    ok &= check(scaled_info.kind == ad::GraphNodeKind::VectorScale, "h * X should create a VectorScale node");
    ok &= check(has_child_kind(scaled_info, ad::GraphNodeKind::ScalarParameter), "VectorScale should keep scalar parameter child");

    ad::Vec sl = X.slice(1, 3);
    ad::GraphInfo slice_info = ad::inspect(sl);
    ok &= check(slice_info.kind == ad::GraphNodeKind::Slice && slice_info.size == 3, "slice should create a Slice node");

    ad::Vec cat = ad::concat(sl, X.slice(0, 1));
    ad::GraphInfo concat_info = ad::inspect(cat);
    ok &= check(concat_info.kind == ad::GraphNodeKind::Concat && concat_info.size == 4, "concat should create a Concat node");

    ad::Vec gathered = ad::gather(X, {0, 4, 2});
    ad::GraphInfo gather_info = ad::inspect(gathered);
    ok &= check(gather_info.kind == ad::GraphNodeKind::Gather && gather_info.size == 3, "gather should create a Gather node");
    ok &= check(ad::count_nodes(gather_info, ad::GraphNodeKind::VectorElement) == 0, "gather should not scalar-lower through Vec::operator[]");

    ad::Vec values = ad::vec_variable("values", 3);
    ad::Vec scattered = ad::scatter_add(values, {0, 4, 4}, 5);
    ad::GraphInfo scatter_info = ad::inspect(scattered);
    ok &= check(scatter_info.kind == ad::GraphNodeKind::ScatterAdd && scatter_info.size == 5, "scatter_add should create a ScatterAdd node");
    ok &= check(ad::count_nodes(scatter_info, ad::GraphNodeKind::VectorElement) == 0, "scatter_add should not scalar-lower through Vec::operator[]");

    ad::Expr s = ad::sum(X);
    ad::GraphInfo sum_info = ad::inspect(s);
    ok &= check(sum_info.kind == ad::GraphNodeKind::Sum, "sum(Vec) should create a Sum reduction node");
    ok &= check(ad::count_nodes(sum_info, ad::GraphNodeKind::VectorElement) == 0, "sum(Vec) should not scalar-lower through Vec::operator[]");

    ad::Expr d = ad::dot(X, X);
    ad::GraphInfo dot_info = ad::inspect(d);
    ok &= check(dot_info.kind == ad::GraphNodeKind::Dot, "dot(Vec, Vec) should create a Dot reduction node");
    ok &= check(ad::count_nodes(dot_info, ad::GraphNodeKind::VectorElement) == 0, "dot(Vec, Vec) should not scalar-lower through Vec::operator[]");

    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 2);
        ad::Vec y = ad::vec_variable("y", 3);
        (void)(x + y);
    }), "mismatched vector binary dimensions should throw");
    ok &= check(throws_cleanly([] {
        ad::DenseMatrix D(2, 2, {1.0, 0.0, 0.0, 1.0});
        ad::Vec x = ad::vec_variable("x", 3);
        (void)(D * x);
    }), "mismatched dense matvec dimensions should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 3);
        (void)x.slice(2, 2);
    }), "invalid slice should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 3);
        (void)ad::gather(x, {0, 3});
    }), "invalid gather index should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 3);
        (void)ad::scatter_add(x, {0, 1}, 4);
    }), "scatter_add index/value size mismatch should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 3);
        (void)ad::scatter_add(x, {0, 1, 4}, 4);
    }), "invalid scatter_add index should throw");

    return ok;
}

bool test_public_headers_do_not_expose_internal_concepts() {
    const std::vector<std::string> headers{
        "ad.h",
        "expr.h",
        "graph_info.h",
        "matrix.h",
        "symbol.h",
        "vec.h",
    };
    const std::vector<std::string> forbidden{
        "builder",
        "plan",
        "emitter",
        "callback",
        "vm",
        "python",
    };

    bool ok = true;
    for (const std::string &header : headers) {
        const std::string body = lower_ascii(read_file(std::string(MOO_AD_SOURCE_DIR) + "/" + header));
        for (const std::string &word : forbidden) {
            ok &= check(body.find(word) == std::string::npos, "forbidden public concept found in " + header + ": " + word);
        }
    }

    const std::string test_cmake = read_file(std::string(MOO_SOURCE_DIR) + "/test/ad/CMakeLists.txt");
    ok &= check(test_cmake.find("test_ad_expr_vec_core") == std::string::npos, "old expr_vec_core test target should not be active");
    ok &= check(test_cmake.find("test_ad_structured_vec") == std::string::npos, "old structured_vec test target should not be active");
    ok &= check(test_cmake.find("test_ad_graph_core") != std::string::npos, "new graph core test should be active");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_scalar_symbol_identity();
    ok &= test_vector_symbol_identity();
    ok &= test_matrix_validation();
    ok &= test_structured_matvec_and_unary_graph();
    ok &= test_slice_concat_scale_and_reductions();
    ok &= test_public_headers_do_not_expose_internal_concepts();
    return ok ? 0 : 1;
}
