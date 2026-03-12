# Quadature-method-numerical-integration
This project compares two numerical integration methods using MATLAB:

- **Gauss Quadrature**
- **Clenshaw-Curtis Quadrature**

The main file evaluates both methods on 6 functions over the interval `[-1,1]`, computes the absolute errors, and plot how the errors change as the number of nodes, n, increases.

**Files
- `Numerical_Integration_Quadrature.m`: main file for running the numerical tests
- `Gaussian_Quadrature.m`: supporting file that implements Gauss quadrature
- `ClenshawCurtis_Quadrature.m`: supporting file that implements Clenshaw-Curtis quadrature

**What this project does
- compare two quadrature methods
- estimate error behaviors for different types of functions
- present numerical convergence through plots and fitted lines

When running this program, open MATLAB, place all three files in the same foler, and run:
**Numerical_Integration_Quadrature.m**
