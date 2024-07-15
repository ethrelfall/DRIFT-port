# DRIFT-port
FEM ports of Fortran DRIFT code (W. Arter, 1994)

The equations are (eqs. 2.2-2.8 in the paper)

$$
n = \phi - \alpha \frac{\partial \phi}{\partial x} - \beta\frac{\partial^3 \phi}{\partial x^3} - \nabla_{\perp}^2 \phi,
$$

$$
\left ( \frac{\partial}{\partial t} + v_E \cdot \nabla_{\perp} \right ) n = - \delta_n \frac{\partial \phi}{\partial x} + \delta_s (y-y_0) \frac{\partial v_{\parallel}}{\partial x} - \nu \nabla_{\parallel}^4 \phi,
$$

$$
\left ( \frac{\partial}{\partial t} + v_E \cdot \nabla_{\perp} \right ) v_{\parallel} = + \delta_s (y-y_0) \frac{\partial \phi}{\partial x} + \left ( \mu_{\perp} \nabla_{\perp}^2 + \mu_{\parallel} \nabla_{\parallel}^2 \right ) v_{\parallel}.
$$

Here

$$
v_E = \left ( -\frac{\partial \phi}{\partial y}, \frac{\partial \phi}{\partial x}, 0 \right ).
$$

Now the game is to reproduce the contour plots in Fig.5.

Boundary conditions: $y$-direction homogeneous Dirichlet, $x$-direction periodic.

Initial data:

Parameter values:

$\nu = \mu_{\perp} = 0.3$

$\mu_{\parallel} \delta_s^2 = 0.03$

$\delta_n = -1$

$\delta_s = 0.2$

$\alpha = 0.744$

$\beta = 0.7$

Box size and centring $L_x =12.5$, $L_y = 50$ and $y_0 = 25$.

Other considerations:
$\nu$ and $\mu_{\perp}$ are started at $0.1$ and cranked up during the simulation.

Equations in weak form:
