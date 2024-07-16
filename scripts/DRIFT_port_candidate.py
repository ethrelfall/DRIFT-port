# DRIFT_port_candidate.py
# Attempt at time-dependent solver of 2D driftwave equations from old DRIFT code
# SU correction (i.e. artificial diffusion) as at
# https://github.com/volpatto/firedrake_scripts/blob/master/scripts/2D/adv-diff_supg_2D.py
# for equations see the repository
# https://github.com/ethrelfall/DRIFT-port/tree/main

from firedrake import *
import math
from irksome import Dt, MeshConstant, TimeStepper, RadauIIA, GaussLegendre, LobattoIIIC
import time
import numpy

meshres = 32
mesh = PeriodicRectangleMesh(meshres, 128, 12.5, 50, direction="x", quadrilateral=True)
V1 = FunctionSpace(mesh, "CG", 1)
V2 = FunctionSpace(mesh, "CG", 1)
V3 = FunctionSpace(mesh, "CG", 1) # space for phi
V3a = FunctionSpace(mesh, "CG", 1)  # space for aux var for phi''
V4 = VectorFunctionSpace(mesh, "CG", 1)  # for drift velocity
V = V1*V2

T = 100  # run for this long
timeres = 1000  # ... in this many steps

t = Constant(0.0)
dt = Constant(T/timeres)

# parameters for irksome
butcher_tableau = RadauIIA(1)

# model parameters
alfa = 0.744
beta = 0.7
delta_n = -1.0
delta_s = 0.2
nu = 0.1
mu_par = 0.75
mu_perp = 0.1
s = 2.0

amp = 1.0  # amplitude of initial data

# numerical parameters
v_eps_sq = 0.00001  # softening parameter for SU correction

x, y = SpatialCoordinate(mesh)
vn = Function(V)
vpar, n = split(vn)
v1, v2 = TestFunctions(V)
driftvel = Function(V4)

# initial data - just give v_parallel an initial profile?
vn.sub(0).interpolate(amp*sin(8*pi*(x-6.25)/12.5)*cos(4*pi*(y-25)/50))
#vn.sub(1).interpolate(amp*sin(8*pi*(x-6.25)/12.5)*cos(5*pi*(y-25)/50))

# compute energy in that initial profile
E = 0.5*assemble( (vn.sub(0)**2)*dx )

print("initial energy=" + str(E))

# TRIALCODE check init data
File("DRIFT_port_v4_candidate.pvd").write(vn.sub(0), vn.sub(1))
#quit()

V3_ext = V3*V3a
phiy3a = TrialFunction(V3_ext)
phi, y3a = split(phiy3a)
phiy3a_s = Function(V3_ext)
phi_s, y3a_s = split(phiy3a_s)
v3, v3a = TestFunctions(V3_ext)

# uses aux variable to handle the third derivative (I did a quick test on that sort of thing for the KdV equation - passed)
Lphi = (inner(grad(phi),grad(v3))+phi*v3-alfa*grad(phi)[0]*v3-beta*grad(y3a)[0]*v3 + (y3a*v3a+inner(grad(phi),grad(v3a))))*dx
Rphi = n*v3*dx
bc1 = DirichletBC(V3_ext.sub(0), 0, [1,2])

h = 0.1*2.0*12.5/meshres  # depends on computational domain size?

# SU damping terms (final two lines) act only on streamline cpt of advection velocity - turned off for now
F = Dt(vpar)*v1*dx + Dt(n)*v2*dx \
  - v1*div(vpar*driftvel)*dx - v2*div(n*driftvel)*dx \
  + (delta_n*grad(phi_s)[0])*(v2)*dx \
  - (delta_s*(y-25.0)*grad(vpar)[0])*(v2)*dx \
  - (delta_s*(y-25.0)*grad(phi_s)[0])*v1*dx \
  + mu_perp*(dot(grad(vpar),grad(v1)))*dx \
  + mu_par*(delta_s**2)*((y-25.0)**2)*(grad(vpar)[0]*(grad(v1)[0]))*dx \
  - nu * inner(grad(phi_s-n-alfa*grad(phi_s)[0]-beta*grad(y3a_s)[0]),grad(v2))*dx \
#  + 0.5*h*(dot(driftvel, grad(vpar)))*dot(driftvel, grad(v1))*(1/((driftvel[0])**2+(driftvel[1])**2+v_eps_sq))*dx \
#  + 0.5*h*(dot(driftvel, grad(n)))*dot(driftvel, grad(v2))*(1/((driftvel[0])**2+(driftvel[1])**2+v_eps_sq))*dx \

# params taken from Cahn-Hilliard example cited above
params = {'snes_monitor': None, 'snes_max_it': 100,
          'snes_linesearch_type': 'l2',
          'ksp_type': 'preonly',
          'pc_type': 'lu', 'mat_type': 'aij',
          'pc_factor_mat_solver_type': 'mumps'}

bc2 = DirichletBC(V.sub(0), 0.0, [1,2])

stepper = TimeStepper(F, butcher_tableau, t, dt, vn, solver_parameters=params, bcs=bc2)

nullspace = VectorSpaceBasis(constant=True)

# this is intended to be direct solver
linparams = {"mat_type": "aij",
          "snes_type": "ksponly",
          "ksp_type": "preonly",
          "pc_type": "lu"}

outfile = File("DRIFT_port_candidate.pvd")

cnt=0
start = time.time()

while float(t) < float(T):
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))
    solve(Lphi==Rphi, phiy3a_s, solver_parameters=linparams, bcs=bc1)
    driftvel.interpolate(as_vector([grad(phi_s)[1],-grad(phi_s)[0]]))
    if(cnt % 20 == 0):
       print("outputting data ...\n")
       vs, ns = vn.split()
       vs.rename("v_par")
       ns.rename("n")
       phiy3a_s.sub(0).rename("phi")
       outfile.write(vs, ns, phiy3a_s.sub(0))
    cnt=cnt+1
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))

end = time.time()
wall_time = end-start

print("done.")
print("\n")
print("wall time:"+str(wall_time)+"\n")
