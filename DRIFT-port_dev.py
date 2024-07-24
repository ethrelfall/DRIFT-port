# DRIFT-port_dev.py
# Attempt at time-dependent solver of 2D driftwave equations from old DRIFT code by Arter
# for which see
# Particle-mesh modelling of drift-wave turbulence, W. Arter, CPC 79 (1994) 381-408
# Attempt made at subtracting k_x=0 part of background from all fields (see calls to "project").
# SU correction (i.e. artificial diffusion) as at
# https://github.com/volpatto/firedrake_scripts/blob/master/scripts/2D/adv-diff_supg_2D.py
# for equations see the repository
# https://github.com/ethrelfall/DRIFT-port/tree/main
# as set up, takes about 8 core-hours (!)

# TODO this seems to need much smaller viscosities nu and mu_perp than Arter's values in order to 
# get a nonzero mode started up - I use 0.01, Arter 0.1 for these params

from firedrake import *
import math
from irksome import Dt, MeshConstant, TimeStepper, RadauIIA, GaussLegendre, LobattoIIIC
import time
import numpy

meshres = 128
mesh = PeriodicRectangleMesh(4*meshres, meshres, 50.0, 12.5, direction="y", quadrilateral=True)
V1 = FunctionSpace(mesh, "CG", 1)
V2 = FunctionSpace(mesh, "CG", 1)
V3 = FunctionSpace(mesh, "CG", 1) # space for phi
V3a = FunctionSpace(mesh, "CG", 1)  # space for aux var for phi'' used to handle third deriv terms
V4 = VectorFunctionSpace(mesh, "CG", 1)  # for drift velocity
V = V1*V2

T=800
timeres=4000

t = Constant(0.0)
dt = Constant(T/timeres)

# parameters for irksome
butcher_tableau = RadauIIA(1)

# model parameters
alfa = 0.744
beta = 0.7
delta_n = -1.0
delta_s = 0.2
nu = 0.01
mu_par = 0.75
mu_perp = 0.01
s = 2.0

amp = 2.0  # amplitude of initial data

# numerical parameters
v_eps_sq = 0.00001  # softening parameter for SU correction

x, y = SpatialCoordinate(mesh)
vn = Function(V)
vpar, n = split(vn)
v1, v2 = TestFunctions(V)
driftvel = Function(V4)

# initial data - just give v_parallel an initial profile
vn.sub(0).interpolate(amp*sin(16*pi*(x-25.0)/50.0)*sin(16*pi*(y-6.25)/12.5)*exp(-(1/0.25**2)*((0.0*(y-6.25)/(0.5*12.5))**2+(((x-0.0)-25)/(0.5*50))**2)))

# TRIALCODE output init data
#File("DRIFT_port_dev.pvd").write(vn.sub(0), vn.sub(1))
#quit()

V3_ext = V3*V3a
phiy3a = TrialFunction(V3_ext)
phi, y3a = split(phiy3a)
phiy3a_s = Function(V3_ext)
phi_s, y3a_s = split(phiy3a_s)
v3, v3a = TestFunctions(V3_ext)

# uses aux variable to handle the third derivative (I did a quick test on that sort of thing for the KdV equation - passed)
Lphi = (inner(grad(phi),grad(v3))+phi*v3-alfa*grad(phi)[1]*v3-beta*grad(y3a)[1]*v3 + (y3a*v3a+grad(phi)[1]*grad(v3a)[1]))*dx
Rphi = n*v3*dx
bc1 = DirichletBC(V3_ext.sub(0), 0, [1,2])

h = 0.1*2.0*12.5/meshres  # depends on computational domain size

# SU damping terms (final two lines) act only on streamline cpt of advection velocity
F = Dt(vpar)*v1*dx + Dt(n)*v2*dx \
  - v1*div(vpar*driftvel)*dx - v2*div(n*driftvel)*dx \
  + (delta_n*grad(phi_s)[1])*(v2)*dx \
  - (delta_s*(x-25.0)*grad(vpar)[1])*(v2)*dx \
  - (delta_s*(x-25.0)*grad(phi_s)[1])*v1*dx \
  + mu_perp*(dot(grad(vpar),grad(v1)))*dx \
  + mu_par*(delta_s**2)*((x-25.0)**2)*(grad(vpar)[1]*(grad(v1)[1]))*dx \
  - nu * inner(grad(phi_s-n-alfa*grad(phi_s)[1]-beta*grad(y3a_s)[1]),grad(v2))*dx \
  + 0.5*h*(dot(driftvel, grad(vpar))+delta_n*grad(phi_s)[0]-delta_s*(x-25.0)*grad(vpar)[1])*dot(driftvel, grad(v1))*(1/sqrt((driftvel[0])**2+(driftvel[1])**2+v_eps_sq))*dx \
  + 0.5*h*(dot(driftvel, grad(n))-delta_s*(x-25.0)*grad(phi_s)[1])*dot(driftvel, grad(v2))*(1/sqrt((driftvel[0])**2+(driftvel[1])**2+v_eps_sq))*dx \

# params taken from Cahn-Hilliard example cited above
params = {'snes_monitor': None, 'snes_max_it': 100,
          'snes_linesearch_type': 'l2',
          'ksp_type': 'preonly',
          'pc_type': 'lu', 'mat_type': 'aij',
          'pc_factor_mat_solver_type': 'mumps'}

bc2 = DirichletBC(V.sub(0), 0.0, [1,2])
bc3 = DirichletBC(V.sub(1), 0.0, [1,2])

stepper = TimeStepper(F, butcher_tableau, t, dt, vn, solver_parameters=params, bcs=bc2)

nullspace = VectorSpaceBasis(constant=True)

# this is intended to be direct solver
linparams = {"mat_type": "aij",
          "snes_type": "ksponly",
          "ksp_type": "preonly",
          "pc_type": "lu"}


# FOR SUBTRACTION
mesh1d = IntervalMesh(4*meshres, 50)
mesh2d = ExtrudedMesh(mesh1d, layers=meshres, layer_height=12.5/meshres, extrusion_type='uniform')
V_sub = FunctionSpace(mesh2d, "CG", 1,  vfamily='R') # for k_x=0 part of phi, to subtract
phi_sub = Function(V_sub)
phi_sub_copyback = Function(V1, name="phi_copyback")
vpar_sub = Function(V_sub)
vpar_sub_copyback = Function(V1)
n_sub = Function(V_sub)
n_sub_copyback = Function(V1)

outfile = File("DRIFT_port_dev.pvd")
energyfile = open("energy.csv", "w")

cnt=0
start = time.time()

while float(t) < float(T):
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))

    #TRIALCODE crank up viscosity at later times?
    #if (float(t) > 60.0):
    #    nu = 0.1
    #    mu_perp = 0.1

    solve(Lphi==Rphi, phiy3a_s, solver_parameters=linparams, bcs=bc1)

    # INTERPOLATE AND DO SUBTRACTION
    phi_sub.interpolate(phi_s)
    U = FunctionSpace(mesh2d, 'CG', 1, vfamily='R', vdegree=0)
    g = Function(U, name='g')
    g.project(phi_sub)

    n_sub.interpolate(vn.sub(1))
    gn = Function(U, name="gn")
    gn.project(n_sub)

    vpar_sub.interpolate(vn.sub(0))
    gvpar = Function(U, name="gvpar")
    gvpar.project(vpar_sub)


    phi_sub_copyback.interpolate(g)
    phiy3a_s.sub(0).interpolate(phi_s-phi_sub_copyback)

    n_sub_copyback.interpolate(gn)
    vn.sub(1).interpolate(vn.sub(1)-n_sub_copyback)

    vpar_sub_copyback.interpolate(gvpar)
    vn.sub(0).interpolate(vn.sub(0)-vpar_sub_copyback)

    driftvel.interpolate(as_vector([grad(phi_s)[1],-grad(phi_s)[0]]))

    if(cnt % 10 == 0):
       print("outputting data ...\n")
       vs, ns = vn.split()
       vs.rename("v_par")
       ns.rename("n")
       phiy3a_s.sub(0).rename("phi")
       outfile.write(vs, ns, phiy3a_s.sub(0), phi_sub_copyback)

    energy = assemble((ns*phiy3a_s.sub(0)+vs*vs)*dx)
    energyfile.write(str(float(t))+", "+str(float(energy))+"\n")

    cnt=cnt+1
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))

end = time.time()
wall_time = end-start

energyfile.close()

print("done.")
print("\n")
print("wall time:"+str(wall_time)+"\n")
