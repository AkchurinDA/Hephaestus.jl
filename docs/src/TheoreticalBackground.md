# Coordinate Transformations

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

# Bisymmetrical Framework Element

## Description

## Elastic Stiffness Matrix

The element elastic stiffness matrix in the local coordinate system is given by:

```math
[\mathbf{k}_{e}^{e}]
= E
\begin{bmatrix}
    +\dfrac{A}{L} & 0 & 0 & 0 & 0 & 0 & 
    -\dfrac{A}{L} & 0 & 0 & 0 & 0 & 0 \\
    0 & +\dfrac{12 I_{zz}}{L^{3}} & 0 & 0 & 0 & +\dfrac{6 I_{zz}}{L^{2}} &
    0 & -\dfrac{12 I_{zz}}{L^{3}} & 0 & 0 & 0 & +\dfrac{6 I_{zz}}{L^{2}} \\
    0 & 0 & +\dfrac{12 I_{yy}}{L^{3}} & 0 & -\dfrac{6 I_{yy}}{L^{2}} & 0 &
    0 & 0 & -\dfrac{12 I_{yy}}{L^{3}} & 0 & -\dfrac{6 I_{yy}}{L^{2}} & 0 \\
    0 & 0 & 0 & +\dfrac{J}{2 (1 + \nu) L} & 0 & 0 & 
    0 & 0 & 0 & -\dfrac{J}{2 (1 + \nu) L} & 0 & 0 \\
    0 & 0 & -\dfrac{6 I_{yy}}{L^{2}} & 0 & +\dfrac{4 I_{yy}}{L} & 0 & 
    0 & 0 & +\dfrac{6 I_{yy}}{L^{2}} & 0 & +\dfrac{2 I_{yy}}{L} & 0 \\
    0 & +\dfrac{6 I_{zz}}{L^{2}} & 0 & 0 & 0 & +\dfrac{4 I_{zz}}{L} & 
    0 & -\dfrac{6 I_{zz}}{L^{2}} & 0 & 0 & 0 & +\dfrac{2 I_{zz}}{L} \\
    -\dfrac{A}{L} & 0 & 0 & 0 & 0 & 0 & 
    +\dfrac{A}{L} & 0 & 0 & 0 & 0 & 0 \\
    0 & -\dfrac{12 I_{zz}}{L^{3}} & 0 & 0 & 0 & -\dfrac{6 I_{zz}}{L^{2}} &
    0 & +\dfrac{12 I_{zz}}{L^{3}} & 0 & 0 & 0 & -\dfrac{6 I_{zz}}{L^{2}} \\
    0 & 0 & -\dfrac{12 I_{yy}}{L} & 0 & +\dfrac{6 I_{yy}}{L^{2}} & 0 &
    0 & 0 & +\dfrac{12 I_{yy}}{L} & 0 & +\dfrac{6 I_{yy}}{L^{2}} & 0 \\
    0 & 0 & 0 & -\dfrac{J}{2 (1 + \nu) L} & 0 & 0 & 
    0 & 0 & 0 & +\dfrac{J}{2 (1 + \nu) L} & 0 & 0 \\
    0 & 0 & -\dfrac{6 I_{yy}}{L^{2}} & 0 & +\dfrac{2 I_{yy}}{L} & 0 & 
    0 & 0 & +\dfrac{6 I_{yy}}{L^{2}} & 0 & +\dfrac{4 I_{yy}}{L} & 0 \\
    0 & +\dfrac{6 I_{zz}}{L^{2}} & 0 & 0 & 0 & +\dfrac{2 I_{zz}}{L} & 
    0 & -\dfrac{6 I_{zz}}{L^{2}} & 0 & 0 & 0 & +\dfrac{4 I_{zz}}{L} \\
\end{bmatrix}
```

The element elastic stiffness matrix in the global coordinate system is given:

```math
[\mathbf{K}_{e}^{e}] = [\mathbf{\Gamma}]^{T} [\mathbf{k}_{e}^{e}] [\mathbf{\Gamma}]
```

The global elastic stiffness matrix is given by:

```math
[\mathbf{K}_{e}] = \text{A}_{i = 1}^{\text{NEL}} [\mathbf{K}_{e}^{i}]
```

## Geometric Stiffness Matrix

The element geometric stiffness matrix in the local coordinate system is given by:

```math
[\mathbf{k}_{g}^{e}]
= \dfrac{F_{x2}}{L}
\begin{bmatrix}
    +1 & 0 & 0 & 0 & 0 & 0 & 
    -1 & 0 & 0 & 0 & 0 & 0 \\
    0 & +\dfrac{6}{5} & 0 & 0 & 0 & +\dfrac{L}{10} &
    0 & -\dfrac{6}{5} & 0 & 0 & 0 & +\dfrac{L}{10} \\
    0 & 0 & +\dfrac{6}{5} & 0 & -\dfrac{L}{10} & 0 &
    0 & 0 & -\dfrac{6}{5} & 0 & -\dfrac{L}{10} & 0 \\
    0 & 0 & 0 & +\dfrac{I_{zz} + I_{yy}}{A} & 0 & 0 &
    0 & 0 & 0 & -\dfrac{I_{zz} + I_{yy}}{A} & 0 & 0 \\
    0 & 0 & -\dfrac{L}{10} & 0 & +\dfrac{2 L^{2}}{15} & 0 &
    0 & 0 & +\dfrac{L}{10} & 0 & -\dfrac{L^{2}}{30} & 0 \\
    0 & +\dfrac{L}{10} & 0 & 0 & 0 & +\dfrac{2 L^{2}}{15} &
    0 & -\dfrac{L}{10} & 0 & 0 & 0 & -\dfrac{L^{2}}{30} \\
    -1 & 0 & 0 & 0 & 0 & 0 & 
    +1 & 0 & 0 & 0 & 0 & 0 \\
    0 & -\dfrac{6}{5} & 0 & 0 & 0 & -\dfrac{L}{10} &
    0 & +\dfrac{6}{5} & 0 & 0 & 0 & -\dfrac{L}{10} \\
    0 & 0 & -\dfrac{6}{5} & 0 & +\dfrac{L}{10} & 0 &
    0 & 0 & +\dfrac{6}{5} & 0 & -\dfrac{L}{10} & 0 \\
    0 & 0 & 0 & -\dfrac{I_{zz} + I_{yy}}{A} & 0 & 0 &
    0 & 0 & 0 & +\dfrac{I_{zz} + I_{yy}}{A} & 0 & 0 \\
    0 & 0 & -\dfrac{L}{10} & 0 & -\dfrac{L^{2}}{30} & 0 &
    0 & 0 & +\dfrac{L}{10} & 0 & +\dfrac{2 L^{2}}{15} & 0 \\
    0 & +\dfrac{L}{10} & 0 & 0 & 0 & -\dfrac{L^{2}}{30} &
    0 & -\dfrac{L}{10} & 0 & 0 & 0 & +\dfrac{2 L^{2}}{15} \\
\end{bmatrix}
```

The element geometric stiffness matrix in the global coordinate system is given:

```math
[\mathbf{K}_{g}^{e}] = [\mathbf{\Gamma}]^{T} [\mathbf{k}_{g}^{e}] [\mathbf{\Gamma}]
```

The global geometric stiffness matrix is given by:

```math
[\mathbf{K}_{g}] = \text{A}_{i = 1}^{\text{NEL}} [\mathbf{K}_{g}^{i}]
```

## Mass Matrix

The element mass matrix in the local coordinate system is given by:

```math
[\mathbf{m}^{e}]
= \dfrac{\rho A L}{420}
\begin{bmatrix}
    +140 & 0 & 0 & 0 & 0 & 0 &
    +70 & 0 & 0 & 0 & 0 & 0 \\
    0 & +156 & 0 & 0 & 0 & +22 L &
    0 & +54 & 0 & 0 & 0 & -13 L \\
    0 & 0 & +156 & 0 & -22 L & 0 &
    0 & 0 & +54 & 0 & +13 L & 0 \\
    0 & 0 & 0 & +\dfrac{140 J}{A} & 0 & 0 &
    0 & 0 & 0 & +\dfrac{70 J}{A} & 0 & 0 \\
    0 & 0 & -22 L & 0 & +4 L^{2} & 0 & 
    0 & 0 & -13 L & 0 & -3 L^{2} & 0 \\
    0 & +22 L & 0 & 0 & 0 & +4 L^{2} &
    0 & +13 L & 0 & 0 & 0 & -3 L^{2} \\ 
    +70 & 0 & 0 & 0 & 0 & 0 & 
    +140 & 0 & 0 & 0 & 0 & 0 \\
    0 & +54 & 0 & 0 & 0 & +13 L & 
    0 & +156 & 0 & 0 & 0 & -22 L \\
    0 & 0 & +54 & 0 & -13 L & 0 & 
    0 & 0 & +156 & 0 & +22 L & 0 \\
    0 & 0 & 0 & +\dfrac{70 J}{A} & 0 & 0 &
    0 & 0 & 0 & +\dfrac{140 J}{A} & 0 & 0 \\
    0 & 0 & +13 L & 0 & -3 L^{2} & 0 & 
    0 & 0 & +22 L & 0 & +4 L^{2} & 0 \\
    0 & -13 L & 0 & 0 & 0 & -3 L^{2} &
    0 & -22 L & 0 & 0 & 0 & +4 L^{2}
\end{bmatrix}
```

The element mass matrix in the global coordinate system is given:

```math
[\mathbf{M}^{e}] = [\mathbf{\Gamma}]^{T} [\mathbf{m}^{e}] [\mathbf{\Gamma}]
```

The global mass matrix is given by:

```math
[\mathbf{M}] = \text{A}_{i = 1}^{\text{NEL}} [\mathbf{M}^{i}]
```

# Analyses Types

## First-Order Elastic Analysis

```math
[\mathbf{K_{e}}] \{\mathbf{U}\} + \{\mathbf{P}\} = \{\mathbf{F}\}
```

## Second-Order Elastic Analysis

```math
[\mathbf{K_{e}} + \mathbf{K_{g}}] \{\mathbf{U}\} + \{\mathbf{P}\} = \{\mathbf{F}\}
```

## First-Order Inelastic Analysis

!!! note
    This feature is currently under development.

## Second-Order Inelastic Analysis

!!! note
    This feature is currently under development.

## Elastic Buckling Analysis

```math
[\mathbf{K_{e}} + \lambda \mathbf{K_{g}}] \{\mathbf{U}\} = \{\mathbf{0}\}
```

## Free Vibration Analysis

```math
[\mathbf{K_{e}} + \omega^{2} \mathbf{M}] \{\mathbf{U}\} = \{\mathbf{0}\}
```