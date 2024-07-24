# KdV_equation.py
# attempt to implement Korteweg-de Vries equation in Firedrake
# motivation: see that higher-order derivatives are handled OK
# compare 1D numerical solution at https://en.wikipedia.org/wiki/Korteweg%E2%80%93De_Vries_equation
# uses auxiliary variables approach for third derivative term, good eyeball agreement with the above solution

from firedrake import *
import math
from irksome import Dt, GaussLegendre, MeshConstant, TimeStepper

meshres = 1000
mesh = PeriodicIntervalMesh(meshres, 2)
V1 = FunctionSpace(mesh, "CG", 1)
V2 = FunctionSpace(mesh, "CG", 1)
V = V1*V2

h=2.0/meshres

T = 20.0
timeres = 2000
t = Constant(0.0)
dt = Constant(T/timeres)

# parameters for irksome
butcher_tableau = GaussLegendre(1)

# model parameters
delta = 0.022

x = SpatialCoordinate(mesh)
solution = Function(V)
u, y = split(solution)
v, z = TestFunctions(V)

# init data
solution.sub(0).interpolate(cos(pi*x[0]))
solution.sub(1).interpolate(-pi*pi*cos(pi*x[0]))

# TRIALCODE output init data
#File("KdV_equation_init.pvd").write(solution.sub(0), solution.sub(1))
#quit()

# final term is SU stabilization of advection term only ...
F = ((Dt(u)*v)*dx) \
   + (u*grad(u)[0]*v)*dx \
   + delta*delta*(grad(y)[0]*v)*dx \
   + (y*z+inner(grad(u),grad(z)))*dx \
   + (0.1*0.5*h*(grad(u)[0])*grad(v)[0])*dx \

# params taken from Cahn-Hilliard example cited above
params = {'snes_monitor': None, 'snes_max_it': 100,
          'snes_linesearch_type': 'l2',
          'ksp_type': 'preonly',
          'pc_type': 'lu', 'mat_type': 'aij',
          'pc_factor_mat_solver_type': 'mumps'}

stepper = TimeStepper(F, butcher_tableau, t, dt, solution, solver_parameters=params)

outfile = File("KdV_equation.pvd")

cnt = 0

while float(t) < float(T):
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))
    if ((cnt%10)==0):
       us, ys = solution.split()
       outfile.write(us)
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))
    cnt=cnt+1

print("done.")
print("\n")
