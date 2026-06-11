// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

bool check(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        return false;
    }
    return true;
}

template <typename Fn>
bool throws_cleanly(Fn &&fn) {
    try {
        fn();
    } catch (const std::exception &) {
        return true;
    }
    return false;
}

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

bool contains_var(const ad::Vars &vars, const ad::Var &var) {
    return std::any_of(vars.values().begin(), vars.values().end(), [&](const ad::Var &candidate) {
        return candidate == var;
    });
}

bool test_map_constructors_and_validation() {
    ad::Vec X = ad::vec_variable("X", 20);

    ad::Map single = ad::Map::single(X, {2, 4});
    ad::Map blocks = ad::Map::blocks(X, 3, 2);
    ad::Map shifted = ad::Map::shifted_blocks(X, 3, 2, 1);
    ad::Map stride = ad::Map::stride(X, 3, 2, 1, 5, 2);
    ad::Map stencil = ad::Map::stencil(X, 3, 0, 5, {0, 4, 7});
    ad::Map explicit_indices = ad::Map::explicit_indices(X, 2, 3, {0, 5, 6, 1, 2, 3});
    ad::Map table = ad::Map::table(X, 2, 3, {7, 4, 0, 8, 9, 10});
    ad::Map expression_source = ad::Map::blocks(X + X, 2, 2);

    bool ok = true;
    ok &= check(single.reps() == 1 && single.local_size() == 2 && single.index(0, 1) == 4, "single map semantics are wrong");
    ok &= check(blocks.reps() == 3 && blocks.local_size() == 2 && blocks.index(2, 1) == 5, "blocks map semantics are wrong");
    ok &= check(shifted.index(2, 1) == 7, "shifted_blocks map semantics are wrong");
    ok &= check(stride.index(2, 1) == 13, "stride map semantics are wrong");
    ok &= check(stencil.index(2, 2) == 17, "stencil map semantics are wrong");
    ok &= check(explicit_indices.index(1, 2) == 3, "explicit_indices map semantics are wrong");
    ok &= check(table.index(1, 0) == 8, "table map semantics are wrong");
    ok &= check(expression_source.source().node_kind() == ad::GraphNodeKind::VectorBinary, "map source should allow structured vector expressions");
    ok &= check(single.kind() == ad::MapKind::Single && single.indices().size() == 2, "single map descriptor metadata is wrong");
    ok &= check(blocks.kind() == ad::MapKind::Blocks && blocks.rep_stride() == 2 && blocks.component_stride() == 1 && blocks.shift() == 0, "blocks map descriptor metadata is wrong");
    ok &= check(shifted.kind() == ad::MapKind::ShiftedBlocks && shifted.rep_stride() == 2 && shifted.shift() == 1, "shifted_blocks map descriptor metadata is wrong");
    ok &= check(stride.kind() == ad::MapKind::Stride && stride.base() == 1 && stride.rep_stride() == 5 && stride.component_stride() == 2, "stride map descriptor metadata is wrong");
    ok &= check(stencil.kind() == ad::MapKind::Stencil && stencil.base() == 0 && stencil.rep_stride() == 5 && stencil.offsets().size() == 3 && stencil.offsets()[2] == 7, "stencil map descriptor metadata is wrong");
    ok &= check(explicit_indices.kind() == ad::MapKind::ExplicitIndices && explicit_indices.indices()[3] == 1, "explicit_indices map descriptor metadata is wrong");
    ok &= check(table.kind() == ad::MapKind::Table && table.indices()[0] == 7, "table map descriptor metadata is wrong");
    ok &= check(!blocks.stores_expanded_indices() && !shifted.stores_expanded_indices() && !stride.stores_expanded_indices() && !stencil.stores_expanded_indices(), "compact map constructors should not store expanded index tables");
    ok &= check(single.stores_expanded_indices() && explicit_indices.stores_expanded_indices() && table.stores_expanded_indices(), "explicit map constructors should report stored index tables");
    ok &= check(throws_cleanly([&] {
        (void)ad::Map::blocks(X, 0, 2);
    }), "zero repetition map should throw");
    ok &= check(throws_cleanly([&] {
        (void)ad::Map::single(X, {20});
    }), "out-of-range map index should throw");
    return ok;
}

bool test_mapped_output_modes_construct() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, x);
    ad::Vec X = ad::vec_variable("X", 6);

    ad::Vec concat = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    });
    ad::Vec sum = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::OutputMode::Sum);
    ad::Vec weighted = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::weighted_sum({1.0, 2.0, 3.0}));
    ad::Vec scatter = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::scatter({0, 2, 2, 3, 1, 3}, 4));

    bool ok = true;
    ok &= check(concat.size() == 6, "concat mapped output size is wrong");
    ok &= check(sum.size() == 2, "sum mapped output size is wrong");
    ok &= check(weighted.size() == 2, "weighted-sum mapped output size is wrong");
    ok &= check(scatter.size() == 4, "scatter mapped output size is wrong");
    ok &= check(throws_cleanly([&] {
        (void)ad::map(local, {ad::bind(x, ad::Map::blocks(X, 3, 2))}, ad::OutputMode::Scatter);
    }), "scatter enum without output indices should throw");
    ok &= check(throws_cleanly([&] {
        (void)ad::map(local, {ad::bind(x, ad::Map::blocks(X, 3, 2))}, ad::MappedOutput::weighted_sum({1.0, 2.0}));
    }), "weighted sum with wrong weight count should throw");
    ok &= check(throws_cleanly([&] {
        (void)ad::map(local, {ad::bind(x, ad::Map::blocks(X, 3, 2))}, ad::MappedOutput::scatter({0, 1}, 2));
    }), "scatter with wrong output index count should throw");
    return ok;
}

bool test_basic_mapped_rhs() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::Vec F = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
        ad::bind(u, ad::Map::blocks(U, 3, 2)),
    });

    ad::GraphInfo info = ad::inspect(F);
    ad::Vars deps = F.variables();
    ad::Function full({X, U}, F);

    bool ok = true;
    ok &= check(F.size() == 6, "mapped RHS concat output size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::MappedFunctionCall, "map should create a MappedFunctionCall node");
    ok &= check(info.children.size() == 2, "mapped inspection should expose source graphs");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorUnary) == 0, "mapped inspection should not inline callee body");
    ok &= check(deps.size() == 12, "mapped RHS variable dependency count is wrong");
    ok &= check(contains_var(deps, X.variables()[0]) && contains_var(deps, U.variables()[5]), "mapped RHS should depend on global sources");
    ok &= check(!contains_var(deps, x.variables()[0]) && !contains_var(deps, u.variables()[0]), "mapped RHS should not leak callee local variables");
    ok &= check(full.output_size() == 6 && full.input_size() == 12, "outer function over mapped RHS has wrong layout");
    return ok;
}

bool test_arbitrary_stencil_local_defect() {
    ad::Vec x = ad::vec_variable("x", 1);
    ad::Vec u = ad::vec_variable("u", 1);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec z = ad::vec_variable("z", 3);
    ad::Vec xc = ad::vec_variable("xc", 1);
    ad::Vec uc = ad::vec_variable("uc", 1);
    ad::Expr h = ad::parameter("h");
    ad::DenseMatrix Dloc(1, 3, {1.0, 2.0, 3.0});
    ad::Function defect({z, xc, uc}, Dloc * z - h * rhs.call({xc, uc}));

    ad::Vec X = ad::vec_variable("X", 20);
    ad::Vec U = ad::vec_variable("U", 4);
    ad::Vec defects = ad::map(defect, {
        ad::bind(z, ad::Map::stencil(X, 4, 0, 4, {0, 4, 7})),
        ad::bind(xc, ad::Map::stride(X, 4, 1, 4, 4, 1)),
        ad::bind(uc, ad::Map::blocks(U, 4, 1)),
    });
    ad::Function full({X, U}, defects);

    bool ok = true;
    ok &= check(defects.size() == 4, "mapped local defect output size is wrong");
    ok &= check(full.parameters().size() == 1 && full.parameters()[0] == h.param(), "mapped local defect should collect callee parameter");
    ok &= check(full.info().output_graph.kind == ad::GraphNodeKind::MappedFunctionCall, "mapped local defect should keep mapped boundary");
    ok &= check(!contains_kind(full.info().output_graph, ad::GraphNodeKind::FunctionCall), "mapped inspection should not inline nested RHS call from callee body");
    ok &= check(contains_var(defects.variables(), X.variables()[19]), "stencil map should include changing local range indices");
    ok &= check(!contains_var(defects.variables(), z.variables()[0]), "mapped local defect should not leak local z variables");
    return ok;
}

bool test_binding_validation_uses_ids() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, x + u);
    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec U = ad::vec_variable("U", 6);

    bool ok = true;
    ok &= check(throws_cleanly([&] {
        (void)ad::map(rhs, {ad::bind(x, ad::Map::blocks(X, 3, 2))});
    }), "missing mapped input binding should throw");
    ok &= check(throws_cleanly([&] {
        (void)ad::map(rhs, {
            ad::bind(x, ad::Map::blocks(X, 3, 2)),
            ad::bind(x, ad::Map::blocks(U, 3, 2)),
        });
    }), "duplicate mapped input binding should throw");
    ok &= check(throws_cleanly([&] {
        ad::Vec x_same_label = ad::vec_variable("x", 2);
        (void)ad::map(rhs, {
            ad::bind(x_same_label, ad::Map::blocks(X, 3, 2)),
            ad::bind(u, ad::Map::blocks(U, 3, 2)),
        });
    }), "same-label different-ID local input binding should throw");
    ok &= check(throws_cleanly([&] {
        (void)ad::map(rhs, {
            ad::bind(x, ad::Map::blocks(X, 3, 2)),
            ad::bind(u, ad::Map::blocks(U, 2, 2)),
        });
    }), "mismatched mapped repetition counts should throw");
    ok &= check(throws_cleanly([&] {
        (void)ad::map(rhs, {
            ad::bind(x, ad::Map::blocks(X, 2, 3)),
            ad::bind(u, ad::Map::blocks(U, 2, 2)),
        });
    }), "mapped local size mismatch should throw");
    return ok;
}

bool test_parameter_sources_and_callee_parameters() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec p_local = ad::vec_parameter("p", 2);
    ad::Function parameterized({x}, x + p_local);

    ad::Vec X = ad::vec_variable("X", 4);
    ad::Vec mapped = ad::map(parameterized, {
        ad::bind(x, ad::Map::blocks(X, 2, 2)),
    });

    ad::Vec P = ad::vec_parameter("P", 4);
    ad::Function source_parameter({x}, ad::sigmoid(x));
    ad::Vec mapped_param_source = ad::map(source_parameter, {
        ad::bind(x, ad::Map::blocks(P, 2, 2)),
    });

    bool ok = true;
    ok &= check(mapped.parameters().size() == 2 && mapped.parameters()[0] == p_local.parameters()[0], "mapped call should collect callee vector parameters");
    ok &= check(mapped_param_source.variables().empty(), "mapped parameter source should not create variable dependencies");
    ok &= check(mapped_param_source.parameters().size() == 4 && mapped_param_source.parameters()[3] == P.parameters()[3], "mapped parameter source should collect selected source parameters");
    return ok;
}

bool test_expression_source_dependencies() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, ad::sigmoid(x));

    ad::Vec X = ad::vec_variable("X", 3);
    ad::DenseMatrix A(4, 3, {1.0, 0.0, 0.0,
                             0.0, 2.0, 0.0,
                             0.0, 0.0, 3.0,
                             1.0, 1.0, 0.0});
    ad::Vec mapped = ad::map(local, {
        ad::bind(x, ad::Map::blocks(A * X, 2, 2)),
    });
    ad::Vars deps = mapped.variables();

    bool ok = true;
    ok &= check(mapped.size() == 4, "mapped expression-source output size is wrong");
    ok &= check(ad::inspect(mapped).children.size() == 1 && ad::inspect(mapped).children[0].kind == ad::GraphNodeKind::DenseMatVec, "mapped expression source should remain structured");
    ok &= check(deps.size() == 3 && contains_var(deps, X.variables()[2]), "mapped expression source should expose source-expression variables");
    ok &= check(!contains_var(deps, x.variables()[0]), "mapped expression source should not leak local input variables");
    return ok;
}

bool test_derivative_and_sparsity_paths_construct() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function rhs({x}, ad::sigmoid(x));
    ad::Vec X = ad::vec_variable("X", 4);
    ad::Vec F = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(X, 2, 2)),
    });
    ad::Vec seed = ad::vec_variable("seed", 4);
    ad::Vec lambda = ad::vec_variable("lambda", F.size());

    bool ok = true;
    ad::Vec jvp = F.forward_diff(ad::vars(X), seed);
    ad::Vec vjp = F.reverse_diff(ad::vars(X), lambda);
    ad::SparsityPattern sp = ad::sparsity(F, ad::vars(X));
    ok &= check(jvp.size() == F.size() && ad::inspect(jvp).kind == ad::GraphNodeKind::MappedFunctionCall, "mapped forward should construct a mapped derivative graph");
    ok &= check(vjp.size() == X.size() && contains_kind(ad::inspect(vjp), ad::GraphNodeKind::MappedFunctionCall), "mapped reverse should construct a mapped derivative graph");
    ok &= check(sp.rows() == F.size() && sp.cols() == X.size(), "mapped sparsity should construct with expected dimensions");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_map_constructors_and_validation();
    ok &= test_mapped_output_modes_construct();
    ok &= test_basic_mapped_rhs();
    ok &= test_arbitrary_stencil_local_defect();
    ok &= test_binding_validation_uses_ids();
    ok &= test_parameter_sources_and_callee_parameters();
    ok &= test_expression_source_dependencies();
    ok &= test_derivative_and_sparsity_paths_construct();
    return ok ? 0 : 1;
}
