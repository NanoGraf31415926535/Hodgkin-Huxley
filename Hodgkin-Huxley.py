import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as animation

# Constants
Cm = 1.0  # membrane capacitance, in uF/cm^2
gNa = 120.0  # maximum conductances, in mS/cm^2
gK = 36.0
gL = 0.3
ENa = 50.0  # Nernst reversal potentials, in mV
EK = -77.0
EL = -54.387

# Define the equations
def alpha_n(V): return 0.01 * (V + 55) / (1 - np.exp(-(V + 55) / 10))
def beta_n(V): return 0.125 * np.exp(-(V + 65) / 80)

def alpha_m(V): return 0.1 * (V + 40) / (1 - np.exp(-(V + 40) / 10))
def beta_m(V): return 4.0 * np.exp(-(V + 65) / 18)

def alpha_h(V): return 0.07 * np.exp(-(V + 65) / 20)
def beta_h(V): return 1 / (1 + np.exp(-(V + 35) / 10))

def dALLdt(X, t):
    V, m, h, n = X
    I_ext = 10 * (t > 10) - 10 * (t > 50)  # External current
    dVdt = (I_ext - gNa * m**3 * h * (V - ENa) - gK * n**4 * (V - EK) - gL * (V - EL)) / Cm
    dmdt = alpha_m(V) * (1 - m) - beta_m(V) * m
    dhdt = alpha_h(V) * (1 - h) - beta_h(V) * h
    dndt = alpha_n(V) * (1 - n) - beta_n(V) * n
    return dVdt, dmdt, dhdt, dndt

# Initial conditions
V0 = -65.0
m0 = alpha_m(V0) / (alpha_m(V0) + beta_m(V0))
h0 = alpha_h(V0) / (alpha_h(V0) + beta_h(V0))
n0 = alpha_n(V0) / (alpha_n(V0) + beta_n(V0))

X0 = [V0, m0, h0, n0]

# Time vector
t = np.linspace(0, 100, 1000)

# Solve ODE
sol = odeint(dALLdt, X0, t)
V = sol[:, 0]
m = sol[:, 1]
h = sol[:, 2]
n = sol[:, 3]

# Additional data for static plots
V_range = np.linspace(-100, 100, 400)

# Equilibrium functions
alpha_n_vals = alpha_n(V_range)
beta_n_vals = beta_n(V_range)
alpha_m_vals = alpha_m(V_range)
beta_m_vals = beta_m(V_range)
alpha_h_vals = alpha_h(V_range)
beta_h_vals = beta_h(V_range)

n_inf = alpha_n_vals / (alpha_n_vals + beta_n_vals)
m_inf = alpha_m_vals / (alpha_m_vals + beta_m_vals)
h_inf = alpha_h_vals / (alpha_h_vals + beta_h_vals)

# Voltage-dependent time constants
tau_n = 1 / (alpha_n_vals + beta_n_vals)
tau_m = 1 / (alpha_m_vals + beta_m_vals)
tau_h = 1 / (alpha_h_vals + beta_h_vals)

# Potassium conductance fit (simulated data)
t_fit = np.linspace(0, 10, 500)
gK_fit = gK * (1 - np.exp(-t_fit / 1))  # Rising phase
gK_fit[t_fit > 5] *= np.exp(-(t_fit[t_fit > 5] - 5) / 1)  # Falling phase

# Stochastic channel activation (simulated data)
np.random.seed(42)
trials = 10
t_stochastic = np.linspace(0, 50, 500)
current_traces = []

for _ in range(trials):
    current_trace = np.random.normal(0, 0.5, len(t_stochastic))  # Simulated random current
    current_trace += np.exp(-t_stochastic / 5)  # Decaying exponential
    current_traces.append(current_trace)

average_current = np.mean(current_traces, axis=0)

# Plotting everything in one figure using subplots
fig, axs = plt.subplots(3, 2, figsize=(15, 12))

# Membrane potential (live plot)
axs[0, 0].set_title('Hodgkin-Huxley Neuron Model')
line, = axs[0, 0].plot(t, V, lw=2)
axs[0, 0].set_xlabel('Time (ms)')
axs[0, 0].set_ylabel('Membrane Potential (mV)')

def update(num, t, V, line):
    line.set_data(t[:num], V[:num])
    return line,

ani = animation.FuncAnimation(fig, update, frames=len(t), fargs=[t, V, line], interval=25, blit=True)

# Equilibrium functions
axs[0, 1].plot(V_range, n_inf, label='$n_{\\infty}$')
axs[0, 1].plot(V_range, m_inf, label='$m_{\\infty}$')
axs[0, 1].plot(V_range, h_inf, label='$h_{\\infty}$')
axs[0, 1].set_xlabel('Voltage (mV)')
axs[0, 1].set_ylabel('Equilibrium value')
axs[0, 1].set_title('Equilibrium functions')
axs[0, 1].legend()

# Voltage-dependent time constants
axs[1, 0].plot(V_range, tau_n, label='$\\tau_n$')
axs[1, 0].plot(V_range, tau_m, label='$\\tau_m$')
axs[1, 0].plot(V_range, tau_h, label='$\\tau_h$')
axs[1, 0].set_xlabel('Voltage (mV)')
axs[1, 0].set_ylabel('Time constant (ms)')
axs[1, 0].set_title('Voltage-dependent time constants')
axs[1, 0].legend()

# Potassium conductance fit
axs[1, 1].plot(t_fit, gK_fit, label='Potassium conductance fit')
axs[1, 1].set_xlabel('Time (ms)')
axs[1, 1].set_ylabel('Conductance (mS/cm^2)')
axs[1, 1].set_title('Potassium conductance fit')
axs[1, 1].legend()

# Stochastic channel activation
for trace in current_traces:
    axs[2, 0].plot(t_stochastic, trace, color='gray', alpha=0.5)
axs[2, 0].plot(t_stochastic, average_current, color='black', label='Average current')
axs[2, 0].set_xlabel('Time (ms)')
axs[2, 0].set_ylabel('Current (nA)')
axs[2, 0].set_title('Stochastic channel activation')
axs[2, 0].legend()

# Hide the empty subplot (2,1)
fig.delaxes(axs[2, 1])

plt.tight_layout()
plt.show()