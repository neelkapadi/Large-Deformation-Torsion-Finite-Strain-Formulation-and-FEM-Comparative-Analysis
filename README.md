# Large-Deformation Torsion: Finite-Strain C++ Formulation and FEM Comparative Analysis

## Overview

This project presents a computational and theoretical investigation of the torsional behaviour of a solid circular shaft undergoing **large rotational deformation**.

A custom **C++ continuum mechanics framework** was developed to evaluate the deformation gradient, strain measures, stress tensors, torsional shear stress, axial normal stress, and torque through numerical integration.

The continuum mechanics results were subsequently validated against a **3D finite element simulation in ANSYS Static Structural** and compared with classical small-deformation torsion theory.

A key objective of the project is to investigate the limitations of classical linear torsion theory under large rotations and demonstrate the resulting **geometric nonlinear effects**, particularly the **Poynting effect**.

---

## Objectives

- Develop a C++-based continuum mechanics formulation for large-deformation torsion.
- Compute analytical and numerical deformation gradients and quantify their error.
- Compare infinitesimal-strain and finite-strain formulations.
- Calculate torsional shear stress and axial normal stress.
- Compute torque through numerical integration over the shaft cross-section.
- Compare finite-strain results with classical torsion theory.
- Validate continuum mechanics predictions using ANSYS FEM.
- Perform a mesh convergence study for the FEM model.

---

## Problem Description

The shaft is subjected to torsion about its longitudinal axis.

The deformation is defined in cylindrical coordinates as:

\[
r = R,\qquad
\theta = \Theta + \alpha Z,\qquad
z = Z
\]

where:

- \(R\) = radial coordinate in the reference configuration
- \(\Theta\) = angular coordinate
- \(Z\) = axial coordinate
- \(\alpha\) = twist per unit length

The corresponding Cartesian deformation map is:

\[
x = X\cos(\alpha Z) - Y\sin(\alpha Z)
\]

\[
y = X\sin(\alpha Z) + Y\cos(\alpha Z)
\]

\[
z = Z
\]

---

## Continuum Mechanics Formulation

### Deformation Gradient

The analytical deformation gradient is:

\[
F =
\begin{bmatrix}
\cos(\alpha Z) & -\sin(\alpha Z) & -\alpha y \\
\sin(\alpha Z) & \cos(\alpha Z) & \alpha x \\
0 & 0 & 1
\end{bmatrix}
\]

The C++ implementation also evaluates the deformation gradient using a **central finite-difference scheme**.

The difference between analytical and numerical results is quantified using the Frobenius norm:

\[
||F-F_{\text{num}}||_F
\]

---

### Strain Measures

Two strain formulations are considered.

#### Infinitesimal Strain

\[
\epsilon =
\frac{1}{2}
\left(
\nabla u + \nabla u^T
\right)
\]

This formulation is used for the small-strain linear elasticity model.

#### Green-Lagrange Strain

\[
E =
\frac{1}{2}(C-I)
\]

where:

\[
C=F^TF
\]

The Green-Lagrange strain is used to account for finite deformation effects.

---

## Material Models

### Small-Strain Linear Elasticity

The Cauchy stress is calculated using:

\[
\sigma^{(s)}
=
\lambda\,\text{tr}(\epsilon)I
+
2\mu\epsilon
\]

where:

\[
\mu =
\frac{E}{2(1+\nu)}
\]

and

\[
\lambda =
\frac{E\nu}
{(1+\nu)(1-2\nu)}
\]

---

### Finite-Strain Elasticity

The second Piola-Kirchhoff stress is calculated as:

\[
S =
\lambda\,\text{tr}(E)I
+
2\mu E
\]

The stress tensor is then transformed into the spatial configuration:

\[
\sigma^{(f)}
=
\frac{1}{J}FSF^T
\]

where:

\[
J=\det(F)
\]

---

## Stress and Torque Calculation

The torsional shear stress is extracted from the Cauchy stress tensor using the circumferential and axial unit vectors:

\[
e_\theta =
(-\sin\theta,\cos\theta,0)
\]

\[
e_z=(0,0,1)
\]

The torsional shear stress is:

\[
\sigma_{z\theta}
=
e_z^T\sigma e_\theta
\]

The axial normal stress is obtained from:

\[
\sigma_{zz}
=
e_z^T\sigma e_z
\]

The total torque is calculated through numerical integration:

\[
T =
\int_A r\sigma_{z\theta}\,dA
\]

The shaft cross-section is discretized into radial and angular material points and the torque is obtained by summing the differential torque contributions.

---

## Classical Torsion Theory

The finite-strain results are compared against the classical small-deformation torsion solution:

\[
\tau(r)=G\alpha r
\]

and

\[
T=GJ_p\alpha
\]

where the polar moment of area for a solid circular shaft is:

\[
J_p=\frac{\pi a^4}{2}
\]

---

## C++ Implementation

The C++ implementation performs the following workflow:

1. Read material, geometry, twist and discretization parameters.
2. Generate material points over the shaft cross-section.
3. Apply the large-rotation deformation map.
4. Calculate analytical deformation gradient.
5. Calculate numerical deformation gradient using central finite differences.
6. Evaluate deformation-gradient error.
7. Calculate Green-Lagrange and infinitesimal strain tensors.
8. Calculate small-strain and finite-strain stress tensors.
9. Extract torsional shear and axial normal stresses.
10. Calculate classical torsion predictions.
11. Numerically integrate torque over the cross-section.
12. Export the complete dataset to CSV for post-processing.

The generated dataset contains deformation, strain, stress, torque contributions, and deformation-gradient error for each discretized material point.

---

## Numerical Parameters

The computational model uses:

| Parameter | Value |
|---|---:|
| Elastic Modulus, \(E\) | 200 GPa |
| Poisson's Ratio, \(\nu\) | 0.3 |
| Shaft Radius, \(a\) | 0.05 m |
| Shaft Length, \(L\) | 0.5 m |
| Twist Rate, \(\alpha\) | 5 rad/m |
| Radial Discretization | 20 |
| Angular Discretization | 36 |
| Finite Difference Step | \(10^{-6}\) |

---

## ANSYS FEM Model

A 3D solid circular shaft was modeled using **ANSYS Static Structural**.

### Boundary Conditions

- One end of the shaft: fully fixed.
- Opposite end: prescribed rotation about the shaft axis.
- Applied end rotation: approximately **2.5 rad (143.24°)**.

### FEM Outputs

The following quantities were extracted:

- Total deformation
- Torsional shear stress
- Axial normal stress
- Reaction torque
- Stress along a radial line

---

## Mesh Convergence Study

Three mesh sizes were investigated:

| Element Size | Nodes | Max Shear Stress (MPa) | Axial Normal Stress (MPa) | Max Deformation (mm) | Max Torque (N·m) |
|---:|---:|---:|---:|---:|---:|
| 10 mm | 22,456 | 19,725 | -562.33 | 94.773 | \(3.7899\times10^6\) |
| 7.5 mm | 40,270 | 19,668 | -551.23 | 94.771 | \(3.7899\times10^6\) |
| 5 mm | 119,189 | 19,606 | -551.87 | 94.768 | \(3.7899\times10^6\) |

The 5 mm mesh was considered converged based on the variation in the extracted quantities.

---

## Results

### Deformation Gradient Validation

The maximum Frobenius-norm error between the analytical and numerical deformation gradients was:
7.18115 x 10^-12

<p align="center">
  <img src="https://github.com/neelkapadi/Large-Deformation-Torsion-Finite-Strain-Formulation-and-FEM-Comparative-Analysis/blob/main/Matlab%20Plots/axial_stress_vs_radius.png" width="800">
</p>

<p align="center">
  <b>Axial Normal Stress vs. Radius</b>
</p>

This demonstrates excellent numerical agreement between the finite-difference implementation and the analytical deformation gradient.

---

### Torque Comparison

For the investigated twist condition:

| Model | Torque |
|---|---:|
| Small-Strain C++ | \(3.77123x10^6\) N·m |
| Finite-Strain C++ | \(4.04570x10^6\) N·m |
| Classical Theory | \(3.77595x10^6\) N·m |
| ANSYS FEM | \(3.78990x10^6\) N·m |

The relative torque error between the finite-strain C++ formulation and ANSYS FEM was:

\[
\boxed{6.32\%}
\]

<p align="center">
  <img src="plots/torque_vs_twist.png" width="800">
</p>

<p align="center">
  <b>Torque vs. Rate of Twist</b>
</p>

<p align="center">
  <img src="plots/shear_stress_vs_radius.png" width="800">
</p>

<p align="center">
  <b>Torsional Shear Stress vs. Radius</b>
</p>

---

## Poynting Effect

The finite-strain formulation predicts non-zero axial normal stresses under large torsional deformation.

This behaviour is associated with the **Poynting effect**, which is not captured by classical small-deformation torsion theory.

The analysis shows that the axial stress varies with radial position and becomes more significant toward the outer surface of the shaft.

This demonstrates the importance of finite-strain kinematics when analyzing shafts subjected to severe rotational deformation.

---

## Key Findings

- The numerical deformation gradient closely matches the analytical solution, with a maximum error of approximately \(7.18\times10^{-12}\).
- Finite-strain theory predicts higher torsional shear stresses than classical small-deformation theory at larger radii.
- Finite-strain formulation captures non-zero axial normal stresses resulting from large rotational deformation.
- The Poynting effect is absent from the classical torsion formulation.
- The finite-strain C++ torque prediction agrees with ANSYS FEM within **6.32% relative error**.
- FEM results show convergence as the mesh is refined from 10 mm to 5 mm.
- Classical relations such as \(\tau=G\alpha r\) and \(T=GJ_p\alpha\) become inadequate for representing the complete stress state under severe rotational deformation.

---
Requirements
C++
C++11 or later
Standard C++ libraries:
<iostream>
<vector>
<cmath>
<fstream>
<iomanip>
<string>
MATLAB

MATLAB is used for post-processing and visualization of the generated CSV data.

ANSYS

ANSYS Static Structural is used for the 3D finite element simulation and mesh convergence study.

---

References
Tadmor, E. B. and Miller, R. E., Modeling Materials: Continuum, Atomistic and Multiscale Techniques, Cambridge University Press.
Sadd, M. H., Continuum Mechanics: Modeling of Material Behavior, Elsevier.
Lai, W. M., Rubin, D., and Krempl, E., Introduction to Continuum Mechanics, 4th ed., Butterworth-Heinemann.

---

Author

Neel M. Kapadi
M.Tech., Department of Mechanical Engineering
Indian Institute of Technology Guwahati

Course: ME 541 – Continuum Mechanics
