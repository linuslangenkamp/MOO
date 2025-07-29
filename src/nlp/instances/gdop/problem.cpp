#include "problem.h"

namespace GDOP {

void FullSweepBuffers::resize(const Mesh& mesh) {
    eval = FixedVector<f64>(mesh.node_count * eval_size);
    jac = FixedVector<f64>(mesh.node_count * jac_size);
    aug_hes = FixedVector<f64>(mesh.node_count * aug_hes_size);
}

void FullSweep::print_jacobian_sparsity_pattern() {
    std::cout << "\n=== LFG Jacobian Sparsity ===\n================================\n";
    for (size_t i = 0; i < lfg.size(); ++i) {
        std::cout << "FunctionLFG[" << i << "] - ";
        if (pc.has_lagrange && i == 0) {
            std::cout << "L - Lagrange term:\n";
        } else if (i >= static_cast<size_t>(pc.f_index_start) && i < static_cast<size_t>(pc.f_index_end)) {
            std::cout << "f[" << (i - pc.f_index_start) << "] - Dynamic Equation:\n";
        } else if (i >= static_cast<size_t>(pc.g_index_start) && i < static_cast<size_t>(pc.g_index_end)) {
            std::cout << "g[" << (i - pc.g_index_start) << "] - Path Constraint:\n";
        } 

        std::cout << "  dx sparsity pattern:\n";
        for (const auto& entry : lfg[i].jac.dx) {
            std::cout << "    x_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
        }

        std::cout << "  du sparsity pattern:\n";
        for (const auto& entry : lfg[i].jac.du) {
            std::cout << "    u_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
        }

        std::cout << "  dp sparsity pattern:\n";
        for (const auto& entry : lfg[i].jac.dp) {
            std::cout << "    p_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
        }

        std::cout << "-----------------------------\n";
    }
    std::cout << "================================\n";
}

void BoundarySweep::print_jacobian_sparsity_pattern() {
    std::cout << "\n=== MR Jacobian Sparsity ===\n================================\n";
    for (size_t i = 0; i < mr.size(); ++i) {
        std::cout << "FunctionMR[" << i << "] - ";
        if (pc.has_mayer && i == 0) {
            std::cout << "M - Mayer term:\n";
        } else if (i >= static_cast<size_t>(pc.r_index_start) && i < static_cast<size_t>(pc.r_index_end)) {
            std::cout << "r[" << (i - pc.r_index_start) << "] - Boundary Constraint:\n";
        } else {
            std::cout << "(unknown):\n";
        }

        std::cout << "  dx0 sparsity pattern:\n";
        for (const auto& entry : mr[i].jac.dx0) {
            std::cout << "    x0_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
        }

        std::cout << "  dxf sparsity pattern:\n";
        for (const auto& entry : mr[i].jac.dxf) {
            std::cout << "    xf_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
        }

        std::cout << "  dp sparsity pattern:\n";
        for (const auto& entry : mr[i].jac.dp) {
            std::cout << "    p_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
        }
        std::cout << "-----------------------------\n";
    }
    std::cout << "================================\n";
}

} // namespace GDOP
