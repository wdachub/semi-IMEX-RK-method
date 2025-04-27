# Runge-Kutta Method Implementation

This repository contains code that implements the Runge-Kutta methods proposed in the paper  
**"Semi-implicit-explicit Runge-Kutta Method for Nonlinear Differential Equations."** 
https://arxiv.org/pdf/2504.09969
It also includes scripts to reproduce several tables from the paper.

## Repository Structure

### `RKmethod/`
Contains implementations of various Runge-Kutta methods, including the semi-implicit-explicit schemes proposed in the paper.

### `test/`
Includes scripts to perform convergence tests for different Runge-Kutta methods, also showing the simplest example of using the function in RKmethod/.

### `tools/`
Contains utility functions used in the tests, such as:
- Construction of differential matrices using the finite difference method
- Conversion of numerical results into LaTeX-formatted tables

### `paper/`
Includes scripts to reproduce Tables 11, 12, and 14 from the paper.

## Citation
If you use this code, please cite the original paper:  
*Semi-implicit-explicit Runge-Kutta Method for Nonlinear Differential Equations.*
```
@article{ding2025semi,
  title={Semi-implicit-explicit Runge-Kutta method for nonlinear differential equations},
  author={Ding, Lingyun},
  journal={arXiv preprint arXiv:2504.09969},
  year={2025}
}
```
