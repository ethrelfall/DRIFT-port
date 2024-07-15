# DRIFT-port
FEM ports of Fortran DRIFT code (W. Arter, 1994)

The equations are

$$
n = \phi - \delta_0 D \phi - \nabla_{\perp}^2 \phi
$$

$$
\left ( \frac{\partial}{\partial t} + v_E \cdot \nabla_{\perp} \right ) n = - \delta_n \frac{\partial \phi}{\partial x} + \delta_s (y-y_0) \frac{\partial v_{\parallel}}{\partial x} - \nu \nabla_{\parallel}^4 \phi
$$

$$
\left ( \frac{\partial}{\partial t} + v_E \cdot \nabla_{\perp} \right ) v_{\parallel} = + \delta_s (y-y_0) \frac{\partial v_{\parallel}}{\partial x} + \left ( \mu_{\perp} \nabla_{\perp}^2 + \mu_{\parallel} \nabla_{\parallel}^2 \right ) v_{\parallel}
$$


