from mpi4py import MPI
from petsc4py import PETSc
import dolfinx
import dolfiny
import basix
import ufl

import gmsh
import numpy as np
import pandas as pd
import time

import gmshapi
import poromechanics.constitutive as pm

save_results = True

"""Mesh from Gmsh"""
name, comm = "drained", MPI.COMM_WORLD
model, gdim, height,elem = gmshapi.mesh_cylindrical_sample()
mesh, subdomains, boundaries = dolfinx.io.gmshio.model_to_mesh(model, comm, rank=0, gdim=3)
omega,gamma_s,gamma_t,gamma_b = 1,2,3,4

"""Function spaces"""
quad_degree, disp_degree = 2, 2

Vue = basix.ufl.element("P", mesh.basix_cell(), disp_degree, shape=(mesh.geometry.dim,))
Vpe = basix.ufl.element("DP", mesh.basix_cell(), 0)
Vqe = basix.ufl.element("RT", mesh.basix_cell(), 1)

Qe = basix.ufl.quadrature_element(mesh.basix_cell(), value_shape=(), degree=quad_degree)
Qeb = basix.ufl.blocked_element(Qe, shape=(mesh.geometry.dim, mesh.geometry.dim), symmetry=True)

Vu = dolfinx.fem.functionspace(mesh, Vue)
Vεp = dolfinx.fem.functionspace(mesh, Qeb)
Vp = dolfinx.fem.functionspace(mesh, Vpe)
Vq = dolfinx.fem.functionspace(mesh, Vqe)

"""Integration measures"""
dx = ufl.Measure("dx", domain=mesh, subdomain_data=subdomains, metadata={"quadrature_degree": quad_degree})
ds = ufl.Measure("ds", domain=mesh, subdomain_data=boundaries, metadata={"quadrature_degree": quad_degree})

"""Solution variables (t=t)"""
u = dolfinx.fem.Function(Vu, name="Displacement")
εp = dolfinx.fem.Function(Vεp, name="Plastic_Strain")
p = dolfinx.fem.Function(Vp, name="Pore_Pressure")
q = dolfinx.fem.Function(Vq, name="Darcy_Flux")

"""Solution variables (t=t-dt)"""
un = dolfinx.fem.Function(Vu)
εpn = dolfinx.fem.Function(Vεp)
pn = dolfinx.fem.Function(Vp)

"""Initial variables (t=0)"""
pi = dolfinx.fem.Function(Vp)

"""Test functions"""
δu = ufl.TestFunction(Vu)
δεp = ufl.TestFunction(Vεp)
δp = ufl.TestFunction(Vp)
δq = ufl.TestFunction(Vq)
δε = ufl.sym(ufl.grad(δu))

"""Constitutive laws"""
functions = {'displacement': u, 'plastic_strain': εp, 'pressure': p,
             'flux': q, 'temperature': None, 'mobilized_friction': None}
old_functions = {'displacement': un, 'plastic_strain': εpn, 'pressure': pn,
                 'temperature': None}
init_conditions = {'pressure': pi, 'temperature': None}

EP = pm.Elastoplasticity(name,mesh,functions,old_functions,init_conditions)

"""Get stress/strain variables"""
ε = EP.total_strain()
εv = EP.volumetric_strain()
εq = EP.elastic_deviator_strain()
εpq = EP.plastic_deviator_strain()
σ = EP.total_stress()
σ_eff = EP.effective_stress()

"""Initialize axial loading (0.1 %/min)"""
load_rate = 1.042e-06 
load_factor = dolfinx.fem.Constant(mesh, 0.0)
plasticity_active = dolfinx.fem.Constant(mesh, 0.0)
strain_active = dolfinx.fem.Constant(mesh, 0.0)

max_strain = 1.5e-2
max_disp = max_strain*height

kk = 25
load = np.linspace(0.064, 1.0, int(kk+1.0))
time_ = load*max_disp/load_rate
dt = time_[1] - time_[0]

"""Variational formulation"""
pc = dolfinx.fem.Constant(mesh,PETSc.ScalarType(-10.0e6))
n = ufl.FacetNormal(mesh)

F0 = ufl.inner(δε,σ)*dx - ufl.inner(pc*n,δu)*ds(gamma_s)
F1 = ufl.inner(δεp, EP.plastic_strain(plasticity_active))*dx
F2 = ufl.inner(δp,EP.mass_balance(dt,strain_active))*dx
F3 = EP.darcy_flux(dt,δq)*dx
F3 += pi*dt*ufl.inner(n,δq)*ds(gamma_t) + pi*dt*ufl.inner(n,δq)*ds(gamma_b)

"""PETSC Nonlinear solver control"""
opts = PETSc.Options(name)
opts["-snes_type"] = "newtonls"
opts["-snes_linesearch_type"] = "basic"
opts["-snes_atol"] = 1.0e-10
opts["-snes_stol"] = 1.0e-6
opts["-snes_rtol"] = 1.0e-12
opts["-snes_max_it"] = 25
opts["-ksp_type"] = "preonly"
opts["-pc_type"] = "lu"
opts["-pc_factor_mat_solver_type"] = "mumps"
opts["-snes_convergence_test"] = "default"
opts["-snes_converged_reason"] = None

"""Coupled problem definition"""
problem = dolfiny.snesblockproblem.SNESBlockProblem([F0,F1,F2,F3], [u,εp,p,q], prefix=name)

"""Essential boundary conditions"""
top_facets = boundaries.indices[boundaries.values == gamma_t]
base_facets = boundaries.indices[boundaries.values == gamma_b]
side_facets = boundaries.indices[boundaries.values == gamma_s]

top_dof_u = dolfinx.fem.locate_dofs_topological(Vu.sub(2),mesh.topology.dim-1,top_facets)
base_dof_u = dolfinx.fem.locate_dofs_topological(Vu.sub(2),mesh.topology.dim-1,base_facets)

side_dof_q = dolfinx.fem.locate_dofs_topological(Vq, mesh.topology.dim-1, side_facets)
q0 = dolfinx.fem.Function(Vq)

"""Results output"""
results = {'σa': [], 'σr': [], 'εa': [], 'εv': [], 'p': [], 't': []}
Volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1.0*dx))
Na = ufl.as_vector([0,0,1])
Nr = ufl.as_vector([np.sqrt(0.5),np.sqrt(0.5),0])

"""Loading loop"""
start = time.time()
for step, factor in enumerate(load):
    """Update axial displacement"""
    load_factor.value = -factor*max_disp
    
    """Update essential bcs"""
    problem.bcs = []
    problem.bcs.append(dolfinx.fem.dirichletbc(load_factor, top_dof_u, Vu.sub(2)))
    problem.bcs.append(dolfinx.fem.dirichletbc(PETSc.ScalarType(0), base_dof_u, Vu.sub(2)))
    problem.bcs.append(dolfinx.fem.dirichletbc(q0,side_dof_q))
    
    if step > 0:
        strain_active.value = 1.0
        plasticity_active.value = 1.0

    print('+++ Step: %0.0f of %0.0f +++' %(step,len(load)-1))
    
    problem.solve()
    assert problem.snes.getConvergedReason() > 0, "Nonlinear solver did not converge!"
    dλ_f = dolfiny.expression.assemble(EP.plastic_consistency(),dx(omega))/Volume
    print('+++ dλ*f = %0.12f +++' %dλ_f)
    
    results['σa'].append(dolfiny.expression.assemble(ufl.dot(σ_eff*Na,Na),dx(omega))/Volume)
    results['σr'].append(dolfiny.expression.assemble(ufl.dot(σ_eff*Nr,Nr),dx(omega))/Volume)
    results['εa'].append(dolfiny.expression.assemble(ufl.dot(ε*Na,Na),dx(omega))/Volume)
    results['εv'].append(dolfiny.expression.assemble(εv,dx(omega))/Volume)
    results['p'].append(dolfiny.expression.assemble(p,dx(omega))/Volume)
    results['t'].append(time_[step])
    
    for source, target in zip([u,εp,p], [un,εpn,pn]):
        with source.vector.localForm() as locs, target.vector.localForm() as loct:
            locs.copy(loct)
            
print("\nSimulation Done! Elapsed time: %0.3f minutes" %((time.time()-start)/60))

for key in results:
        results[key].insert(0,0)

σa = np.array(results['σa'])/-1e6
σr = np.array(results['σr'])/-1e6
εa = np.array(results['εa'])*-1e2
εv = np.array(results['εv'])*-1e2
p = np.array(results['p'])/1e6
t = np.array(results['t'])/60

if save_results:
    dict_ = {'Axial_stress':σa,'Radial_stress':σr,'Axial_strain':εa,'Volumetric_strain':εv,'Pore_pressure':p,'time':t}
    df = pd.DataFrame(dict_,dtype=np.float64)
    df.to_csv(f"{name}.csv")