# --- simulation of equivalent circuit: "[C1-P2,R3]"
using PyCall, RobustModels, EquivalentCircuits
circuits = PyCall.pyimport("impedance.models.circuits")

# --- Python Functions -----------------------------------------------------------------------------------------------------
# --- Phyton code ----------------------------------------------------------------------------------------------------------
py"""
# -------------------------------------------------------------------------
# origin: Michael Grubmüller:
# simulation of equivalent circuit: "[C1-P2,R3]"
# -------------------------------------------------------------------------
import numpy as np

def zR(f, R):
    zR = R * np.ones(len(f))
    return zR

def zC(f, C):
    zC = 1/(1j*2*np.pi*f*C)
    return zC

def zCPE(f, Q, alpha):
    zCPE = 1/((1j*2*np.pi*f)**alpha * Q)
    return zCPE

    # parallel elements
def zpar(z1, z2):
    zpar = z1*z2 / (z1+z2)
    return zpar

def SimulateEQ(f, C1, CPE1_Q, CPE1_alpha, R1):
    Z_sim = zpar(zC(f,C1) + zCPE(f,CPE1_Q, CPE1_alpha), zR(f,R1))
    return Z_sim

# Frequency range
f = np.logspace(np.log10(0.1), np.log10(100e3))

# Component parameters
C1          = 6.76e-1
CPE1_Q      = 3.75
CPE1_alpha  = 0.24
R1          = 50.018

# Resulting impedance
Z_simulated = zpar(zC(f,C1) + zCPE(f,CPE1_Q, CPE1_alpha), zR(f,R1))

""" # --- end Python-Block -------------------------------------------------------------------------------------------------

# --- Julia functions: -----------------------------------------------------------------------------------------------------
function Z_R(_f::Vector{Float64}, _R::Number)
    return _R .* ones(size(_f))
end

function Z_C(_f::Vector{Float64}, _C::Number)
    return 1.0 ./ (1im * 2 * pi * _C .* _f)
end

function Z_P(_f::Vector{Float64}, _Pw::Number, _Pn::Number)
    return 1.0 ./ (_Pw .* (1im * 2 * pi .* _f).^_Pn)
end

function parECE(_z1, _z2)
    return _z1 .* _z2 ./ (_z1 .+ _z2)
end

# --- variables / parameters: -------------------------------------------------------------------------
n_elements = 10; f1 = 1e-1; fend = 1e+5;
f_      = exp10.(LinRange(log10(f1), log10(fend), n_elements))
C1      = 6.76e-1
P2w     = 3.75
P2n     = 0.24
R3      = 50.018
circuitStr_ImpPy = "p(C1-CPE2,R3)"
param_ImpPy      = [C1, P2w, P2n, R3]
circuitStr_ECircJL = "[C1-P2,R3]"
param_ECircJL = (C1 = C1, P2w = P2w, P2n = P2n, R3 = R3)

# --- Julia Calculation Block ----------------------------------------------------------------------------------------------
zR3_ = Z_R(f_, R3)
zC1_ = Z_C(f_, C1)
z_P2_ = Z_P(f_, P2w, P2n)

# --- simulation:
Z_sim_own = parECE(zC1_ + z_P2_, zR3_)

# --- impedance.py Calculation Block ---------------------------------------------------------------------------------------
equiv_circuit_PyObj = circuits.CustomCircuit(initial_guess= param_ImpPy, circuit= circuitStr_ImpPy)
Z_simulated_ImpPy = equiv_circuit_PyObj.predict(f_, use_initial = true)

# --- call own python function ---------------------------------------------------------------------------------------------
Z_Py = py"SimulateEQ"(f_, C1, P2w, P2n, R3)

# --- simulation via package EquivalentCircuits-Package: -------------------------------------------------------------------
circfunc_ = EquivalentCircuits.circuitfunction(circuitStr_ECircJL)
Z_ECircJL = EquivalentCircuits.simulateimpedance_noiseless(circfunc_, param_ECircJL, f_)

# --- comparison of results: -----------------------------------------------------------------------------------------------
ΔQz_ImpPy_Qz_own     = RobustModels.mean(abs.(Z_simulated_ImpPy - Z_sim_own))
ΔQz_ImpPy_Qz_Py      = RobustModels.mean(abs.(Z_simulated_ImpPy - Z_Py))
ΔQz_ImpPy_Qz_ECircJL = RobustModels.mean(abs.(Z_simulated_ImpPy - Z_ECircJL))
