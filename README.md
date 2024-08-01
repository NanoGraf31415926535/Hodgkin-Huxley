Here's a README file for the Hodgkin-Huxley model simulation. This file will provide an overview of the project, instructions for setting it up, and explanations for the various plots.

### README.md
```markdown
# Hodgkin-Huxley Model Simulation

This project simulates the Hodgkin-Huxley model of neuronal activity using Python. The model is visualized through several plots, including a live-updating membrane potential plot, and static plots of equilibrium functions, voltage-dependent time constants, potassium conductance fit, and stochastic channel activation.

## Description

The Hodgkin-Huxley model describes how action potentials in neurons are initiated and propagated. This simulation solves the model's differential equations and provides a dynamic visualization of the results.

## Features

- **Live-updating membrane potential plot**: Shows how the membrane potential changes over time.
- **Equilibrium functions plot**: Displays the equilibrium values of gating variables \(m\), \(n\), and \(h\) as functions of voltage.
- **Voltage-dependent time constants plot**: Shows the time constants for the gating variables as functions of voltage.
- **Potassium conductance fit plot**: Fits a model to the potassium conductance data.
- **Stochastic channel activation plot**: Simulates and displays stochastic activation of ion channels.

## Usage

To run the simulation, execute the following command:
```bash
python hodgkin_huxley_simulation.py
```

The script will open a window displaying the following plots:

1. **Live-updating membrane potential plot**:
   - Shows the change in membrane potential over time as the simulation progresses.
   
2. **Equilibrium functions plot**:
   - Displays \(m_{\infty}(V)\), \(n_{\infty}(V)\), and \(h_{\infty}(V)\) as functions of voltage.
   
3. **Voltage-dependent time constants plot**:
   - Shows the time constants \( \tau_m(V) \), \( \tau_n(V) \), and \( \tau_h(V) \) as functions of voltage.
   
4. **Potassium conductance fit plot**:
   - Fits and shows the potassium conductance \( g_K n^4(t) \) based on the model.

5. **Stochastic channel activation plot**:
   - Simulates and plots the stochastic activation of ion channels over multiple trials and their average.

## Code Explanation

### Differential Equations

The core of the simulation solves the Hodgkin-Huxley equations:
```python
def dALLdt(X, t):
    V, m, h, n = X
    I_ext = 10 * (t > 10) - 10 * (t > 50)
    dVdt = (I_ext - gNa * m**3 * h * (V - ENa) - gK * n**4 * (V - EK) - gL * (V - EL)) / Cm
    dmdt = alpha_m(V) * (1 - m) - beta_m(V) * m
    dhdt = alpha_h(V) * (1 - h) - beta_h(V) * h
    dndt = alpha_n(V) * (1 - n) - beta_n(V) * n
    return dVdt, dmdt, dhdt, dndt
```

### Plots

- **Live-updating plot**: Utilizes `matplotlib.animation.FuncAnimation` to update the membrane potential in real-time.
- **Static plots**: Generated using `matplotlib` to display equilibrium functions, time constants, potassium conductance fit, and stochastic channel activation.

## References

- Hodgkin, A. L., & Huxley, A. F. (1952). A quantitative description of membrane current and its application to conduction and excitation in nerve. *The Journal of Physiology*, 117(4), 500-544.
- Patlak, J. B., & Ortiz, M. (1985). Two modes of gating during late Na+ channel currents in frog sartorius muscle. *The Journal of General Physiology*, 85(1), 31-51.
https://neuronaldynamics.epfl.ch/online/Ch2.S2.html
https://en.wikipedia.org/wiki/Hodgkin–Huxley_model

