- ``\rho`` - rotation about the local ``y``-axis, bringing the ``x_{1}``-axis into the ``x_{1}``-``y`` plane. This angle is computed automatically.
- ``\chi`` - rotation about the ``z_{1}``-axis, bringing the ``x_{1}``-axis into coincidence with the global ``x^{\prime}``-axis. This angle is computed automatically.
- ``\omega`` - rotation about the ``x_{2}``-axis, bringing the ``y_{2}``-axis into coincidence with the global ``y^{\prime}``-axis and the ``z_{2}``-axis into coincidence with the global ``z^{\prime}``-axis. This angle must be provided by the user.

```math
\begin{align*}
    \{\mathbf{x}^{\prime}\}
    &= [\mathbf{\gamma}] \{\mathbf{x}\} \\
    &= [\mathbf{\gamma}_{\omega}] [\mathbf{\gamma}_{\chi}] [\mathbf{\gamma}_{\rho}] \{\mathbf{x}\} \\
    &= 
    \begin{bmatrix}
        1 & 0 & 0 \\
        0 & +\cos{\omega} & +\sin{\omega} \\
        0 & -\sin{\omega} & +\cos{\omega}
    \end{bmatrix}
    \begin{bmatrix}
        +\cos{\chi} & +\sin{\chi} & 0 \\
        -\sin{\chi} & +\cos{\chi} & 0 \\
        0 & 0 & 1
    \end{bmatrix}
    \begin{bmatrix}
        +\cos{\rho} & 0 & +\sin{\rho} \\
        0 & 1 & 0 \\
        -\sin{\rho} & 0 & +\cos{\rho}
    \end{bmatrix}
    \{\mathbf{x}\}
\end{align*}
```

```math
[\mathbf{\Gamma}] = 
\begin{bmatrix}
    [\mathbf{\gamma}] & [\mathbf{0}] & [\mathbf{0}] & [\mathbf{0}] \\
    [\mathbf{0}] & [\mathbf{\gamma}] & [\mathbf{0}] & [\mathbf{0}] \\
    [\mathbf{0}] & [\mathbf{0}] & [\mathbf{\gamma}] & [\mathbf{0}] \\
    [\mathbf{0}] & [\mathbf{0}] & [\mathbf{0}] & [\mathbf{\gamma}]
\end{bmatrix}
```

