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
\left ( \frac{\partial}{\partial t} + v_E \cdot \nabla_{\perp} \right ) v_{\parallel} = + \delta_s (y-y_0) \frac{\partial v_{\parallel}}{\partial x} + \left ( \mu_{\perp} \nabla_{\perp}^2 + \mu_{\parallel} \nabla_{\parallel}^2 \right ) v_{\parallel}.
$$

Here

$$
v_E = \left ( -\frac{\partial \phi}{\partial y}, \frac{\partial \phi}{\partial x}, 0 \right ).
$$

Now the game is to reproduce the contour plots in Fig.5.

Initial data:
Parameter values:
Other considerations:
