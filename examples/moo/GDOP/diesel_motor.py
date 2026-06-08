# SPDX-License-Identifier: LGPL-3.0-or-later
#
# This file is part of MOO - Modelica / Model Optimizer
# Copyright (C) 2026 University of Applied Sciences and Arts
# Bielefeld, Faculty of Engineering and Mathematics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from math import pi
from pathlib import Path

from moo import gdop_model


def sqrt(x):
    return x ** 0.5


p_amb = 1.0111134146341463e005
T_amb = 2.9846362195121958e002
gamma_a = 1.3964088397790055e000
R_a = 2.8700000000000000e002
cp_a = 1.0110000000000000e003
cv_a = 7.2400000000000000e002
p_es = 1.0111134146341463e005
gamma_e = 1.2734225621414914e000
R_e = 2.8600000000000000e002
cp_e = 1.3320000000000000e003
Hlhv = 4.2900000000000000e007
gamma_cyl = 1.3500000000000001e000
T_im = 3.0061857317073162e002
n_cyl = 6
V_D = 1.2699999999999999e-002
r_c = 1.7300000000000001e001
A_wg_eff = 8.8357293382212933e-004
state_norm1 = 2.2000000000000000e002
state_norm2 = 2.0000000000000000e005
state_norm3 = 3.0000000000000000e005
state_norm4 = 1.0000000000000000e004
control_norm1 = 2.5000000000000000e002
control_norm3 = 2.5400000000000000e005
Psi_max = 1.4374756793329366e000
dot_m_c_corr_max = 5.2690636559024850e-001
eta_igch = 6.8768988621327665e-001
c_fr1 = 7.1957840228405290e-001
c_fr2 = -1.4144357053459333e-001
c_fr3 = 3.5904197283929118e-001
eta_sc = 1.0515746242284574e000
x_cv = 5.7966369236798054e-001
A_t_eff = 9.9877716035454514e-004
eta_c = 5.2227332577901808e-001
eta_t = 6.8522930965034778e-001
eta_vol = 8.9059680994120261e-001
w_fric = 2.4723010996875069e-005


def diesel_terms(w_ice, p_im, p_em, w_tc, u_f, u_wg):
    W_ICE = state_norm1 * w_ice
    P_IM = state_norm2 * p_im
    P_EM = state_norm3 * p_em
    W_TC = state_norm4 * w_tc
    U_F = control_norm1 * u_f

    Pi_c = P_IM / p_amb
    w_tc_corr = state_norm4 * w_tc
    Pi_c_max = ((((w_tc_corr ** 2) * (0.04 ** 2) * Psi_max) / (2 * cp_a * T_amb)) + 1) ** (gamma_a / (gamma_a - 1))
    Cm_temp = 1 - ((Pi_c / Pi_c_max) ** 2)
    dot_m_c = dot_m_c_corr_max * sqrt(Cm_temp)
    P_c = dot_m_c * cp_a * T_amb * ((Pi_c ** ((gamma_a - 1) / gamma_a)) - 1) / eta_c

    dot_m_ci = eta_vol * P_IM * W_ICE * V_D / (4 * pi * R_a * T_im)
    dot_m_f = U_F * W_ICE * n_cyl * 1e-6 / (4 * pi)

    eta_ig = eta_igch * (1 - (1 / (r_c ** (gamma_cyl - 1))))
    T_pump = V_D * (P_EM - P_IM)
    T_ig = n_cyl * Hlhv * eta_ig * u_f * control_norm1 * 1e-6
    speed_krpm = W_ICE * 60 / (2 * pi * 1000)
    T_fric = V_D * 1e5 * (c_fr1 * (speed_krpm ** 2) + c_fr2 * speed_krpm + c_fr3)
    T_ice = (T_ig - T_fric - T_pump) / (4 * pi)

    Pi_e = P_EM / P_IM
    q_in = dot_m_f * Hlhv / (dot_m_f + dot_m_ci)
    x_p = 1 + q_in * x_cv / (cv_a * T_im * (r_c ** (gamma_a - 1)))
    T_eo = (
        eta_sc
        * (Pi_e ** (1 - 1 / gamma_a))
        * (r_c ** (1 - gamma_a))
        * (x_p ** (1 / gamma_a - 1))
        * (q_in * ((1 - x_cv) / cp_a + x_cv / cv_a) + T_im * (r_c ** (gamma_a - 1)))
    )

    Pi_t = p_es / P_EM
    Pi_ts = sqrt(Pi_t)
    Psi_t = sqrt((2 * gamma_e / (gamma_e - 1)) * ((Pi_ts ** (2 / gamma_e)) - (Pi_ts ** ((gamma_e + 1) / gamma_e))))
    dot_m_t = P_EM * Psi_t * A_t_eff / sqrt(T_eo * R_e)
    P_t = dot_m_t * cp_e * T_eo * eta_t * (1 - (Pi_t ** ((gamma_e - 1) / gamma_e)))

    Pi_wg = p_amb / P_EM
    Psi_wg = sqrt(2 * gamma_e / (gamma_e - 1) * ((Pi_wg ** (2 / gamma_e)) - (Pi_wg ** ((gamma_e + 1) / gamma_e))))
    dot_m_wg = P_EM * Psi_wg * A_wg_eff * u_wg / sqrt(T_eo * R_e)

    return {
        "T_ice": T_ice,
        "Cm_temp": Cm_temp,
        "dot_m_ci": dot_m_ci,
        "T_eo": T_eo,
        "dot_m_f": dot_m_f,
        "dot_m_t": dot_m_t,
        "dot_m_wg": dot_m_wg,
        "P_t": P_t,
        "P_c": P_c,
        "W_TC": W_TC,
    }


def rhs(w_ice, p_im, p_em, w_tc, u_f, u_wg):
    t = diesel_terms(w_ice, p_im, p_em, w_tc, u_f, u_wg)
    return [
        0.0012987012987013 * t["T_ice"],
        20.2505119361145 * ((0.526906365590249 * sqrt(t["Cm_temp"])) - t["dot_m_ci"]),
        0.0476078551344513 * (t["T_eo"] * (t["dot_m_ci"] + t["dot_m_f"] - t["dot_m_t"] - t["dot_m_wg"])),
        0.0001 * ((t["P_t"] - t["P_c"]) / (0.000197779559297041 * t["W_TC"]) - w_fric * t["W_TC"] * t["W_TC"]),
    ]


def mayer(w_ice, p_im, p_em, w_tc):
    return (
        (w_ice - 0.515309170685596) ** 2
        + (p_im - 0.547055854225991) ** 2
        + (p_em - 0.381048005791294) ** 2
        + (w_tc - 0.271443000537680) ** 2
    )


def build_model(name: str = "DieselMotor"):
    model = gdop_model(name)

    w_ice = model.add_state("w_ice", start=2.4989941562646081e-01, lb=4 / state_norm1, ub=220 / state_norm1)
    p_im = model.add_state("p_im", start=5.0614999999999999e-01, lb=0.8 * p_amb / state_norm2, ub=2 * p_amb / state_norm2)
    p_em = model.add_state("p_em", start=3.3926666666666666e-01, lb=p_amb / state_norm3, ub=3 * p_amb / state_norm3)
    w_tc = model.add_state("w_tc", start=6.8099999999999994e-02, lb=300 / state_norm4, ub=10000 / state_norm4)
    u_f = model.add_control("u_f", lb=0.0, ub=250 / control_norm1, guess=0.5)
    u_wg = model.add_control("u_wg", lb=0.0, ub=1.0, guess=0.5)

    dynamics = rhs(w_ice, p_im, p_em, w_tc, u_f, u_wg)
    for state, state_rhs in zip([w_ice, p_im, p_em, w_tc], dynamics):
        model.set_dynamics(state, state_rhs)

    model.add_mayer(mayer(w_ice, p_im, p_em, w_tc))
    model.add_lagrange(diesel_terms(w_ice, p_im, p_em, w_tc, u_f, u_wg)["dot_m_f"])
    model.set_time_fixed(t0=0.0, tf=0.5)
    model.mesh(intervals=50, nodes=3)
    model.mesh_refinement(2, 5)
    model.solver(tolerance=1e-12)
    return model


if __name__ == "__main__":
    out = Path("build/moo/GDOP/diesel_motor")
    run = build_model().run(out)
    run.result.plot.all(save=out / "solution.png")
    raise SystemExit(run.returncode)
