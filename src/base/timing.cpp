// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2025 University of Applied Sciences and Arts
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

#include <base/timing.h>
#include <base/util.h>
#include <base/log.h>

TimingNode::TimingNode(std::string n, TimingNode* p)
    : name(std::move(n)), parent(p) {}

TimingNode* TimingNode::add_child(std::unique_ptr<TimingNode> child) {
    children.push_back(std::move(child));
    children.back()->parent = this;
    return children.back().get();
}

void TimingNode::finalize() {};

TimingTree& TimingTree::instance() {
    static TimingTree inst;
    return inst;
}

TimingNode* TimingTree::current_node() {
    return current;
}

void TimingTree::set_current(TimingNode* n) {
    current = n;
}

TimingNode* TimingTree::add_child(std::unique_ptr<TimingNode> child) {
    child->parent = current;
    current->children.push_back(std::move(child));
    return current->children.back().get();
}

TimingTree::TimingTree() {
    root = std::make_unique<TimingNode>("root");
    current = root.get();
}

CountedTimingNode::CountedTimingNode(std::string n, TimingNode* p, int cnt)
    : TimingNode(n, p), cnt(cnt) {}


// ==== printing ====

void TimingNode::print_table_tree(const FixedTableFormat<4>& ftf, const std::string& prefix, bool is_last) const {
    f64 duration_ms = nano_to_ms(duration_nano);

    std::string count_str = "-";
    std::string avg_str = "-";

    if (const auto* counted = dynamic_cast<const CountedTimingNode*>(this)) {
        count_str = fmt::format("{}", counted->cnt);
        if (counted->cnt > 0) {
            avg_str = fmt::format("{:.3f}", duration_ms / counted->cnt);
        }
    }

    std::string branch = " ";
    std::string duration_string = "-";
    if (parent != nullptr) {
        branch = is_last ? "└─ " : "├─ ";
        duration_string = fmt::format("{:.3f}", duration_ms);
    }

    std::string tree_col = prefix + branch + name;

    const size_t tree_col_width = ftf.col_widths[0];
    if (tree_col.length() > tree_col_width) {
        tree_col = tree_col.substr(0, tree_col_width - 3) + "...";
    }
    else {
        tree_col += std::string(tree_col_width - tree_col.length(), ' ');
    }

    Log::row(ftf, tree_col, duration_string, count_str, avg_str);

    for (size_t ch_idx = 0; ch_idx < children.size(); ch_idx++) {
        children[ch_idx]->print_table_tree(ftf, prefix + (is_last ? "   " : "│  "), ch_idx == children.size() - 1);
    }
}

void TimingTree::print_tree_table() const {
    if (!root) return;

    FixedTableFormat<4> ftf({60, 12, 8, 12}, {Align::Left, Align::Right, Align::Right, Align::Right});

    Log::start_module(ftf, "Walltime Overview");
    Log::row(ftf, std::string(27, ' ') + "Section" + std::string(26, ' '), "Total [ms]", "Count", "Avg [ms]");
    Log::dashes(ftf);
    root->print_table_tree(ftf);
    Log::dashes_ln(ftf);
}
