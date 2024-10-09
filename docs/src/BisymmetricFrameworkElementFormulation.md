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