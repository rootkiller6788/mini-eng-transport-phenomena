#include "cfd_types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ===== Grid 1D ===== */
CfdGrid1D *cfd_grid1d_create(int nx, double L) {
    if (nx<2||L<=0.0) return NULL;
    CfdGrid1D*g=malloc(sizeof(CfdGrid1D));if(!g)return NULL;
    g->nx=nx;g->nnodes=nx+1;g->L=L;g->dx=L/(double)nx;
    g->x=malloc((size_t)g->nnodes*sizeof(double));
    g->xc=malloc((size_t)nx*sizeof(double));
    if(!g->x||!g->xc){free(g->x);free(g->xc);free(g);return NULL;}
    for(int i=0;i<g->nnodes;i++)g->x[i]=(double)i*g->dx;
    for(int i=0;i<nx;i++)g->xc[i]=((double)i+0.5)*g->dx;
    return g;
}
void cfd_grid1d_destroy(CfdGrid1D*g){if(g){free(g->x);free(g->xc);free(g);}}

/* ===== Grid 2D ===== */
CfdGrid2D *cfd_grid2d_create(int nx,int ny,double Lx,double Ly){
    if(nx<2||ny<2||Lx<=0.0||Ly<=0.0)return NULL;
    CfdGrid2D*g=malloc(sizeof(CfdGrid2D));if(!g)return NULL;
    g->nx=nx;g->ny=ny;g->ncells=nx*ny;g->nnodes=(nx+1)*(ny+1);
    g->Lx=Lx;g->Ly=Ly;g->dx=Lx/(double)nx;g->dy=Ly/(double)ny;
    g->topology=CFD_TOPOLOGY_STRUCTURED;g->cell_shape=CFD_CELL_QUAD;
    g->x=g->y=g->xc=g->yc=g->cell_vol=g->face_area_x=g->face_area_y=NULL;
    return g;
}
int cfd_grid2d_allocate_all(CfdGrid2D*g){
    if(!g)return -1;
    int nc=g->ncells,nn=g->nnodes,nfx=(g->nx+1)*g->ny,nfy=g->nx*(g->ny+1);
    g->x=malloc((size_t)nn*sizeof(double));g->y=malloc((size_t)nn*sizeof(double));
    g->xc=malloc((size_t)nc*sizeof(double));g->yc=malloc((size_t)nc*sizeof(double));
    g->cell_vol=malloc((size_t)nc*sizeof(double));
    g->face_area_x=malloc((size_t)nfx*sizeof(double));
    g->face_area_y=malloc((size_t)nfy*sizeof(double));
    if(!g->x||!g->y||!g->xc||!g->yc||!g->cell_vol||!g->face_area_x||!g->face_area_y)
        {cfd_grid2d_destroy(g);return -1;}
    double dx=g->dx,dy=g->dy;
    for(int j=0;j<=g->ny;j++)for(int i=0;i<=g->nx;i++)
        {int idx=CFD_IDX2(i,j,g->nx+1);g->x[idx]=(double)i*dx;g->y[idx]=(double)j*dy;}
    for(int j=0;j<g->ny;j++)for(int i=0;i<g->nx;i++)
        {int idx=CFD_IDX2(i,j,g->nx);g->xc[idx]=((double)i+0.5)*dx;g->yc[idx]=((double)j+0.5)*dy;g->cell_vol[idx]=dx*dy;}
    for(int i=0;i<nfx;i++)g->face_area_x[i]=dy;
    for(int i=0;i<nfy;i++)g->face_area_y[i]=dx;
    return 0;
}
void cfd_grid2d_destroy(CfdGrid2D*g){if(g){free(g->x);free(g->y);free(g->xc);free(g->yc);free(g->cell_vol);free(g->face_area_x);free(g->face_area_y);free(g);}}

/* ===== Grid 3D ===== */
CfdGrid3D *cfd_grid3d_create(int nx,int ny,int nz,double Lx,double Ly,double Lz){
    if(nx<2||ny<2||nz<2||Lx<=0.0||Ly<=0.0||Lz<=0.0)return NULL;
    CfdGrid3D*g=malloc(sizeof(CfdGrid3D));if(!g)return NULL;
    g->nx=nx;g->ny=ny;g->nz=nz;g->ncells=nx*ny*nz;g->nnodes=(nx+1)*(ny+1)*(nz+1);
    g->Lx=Lx;g->Ly=Ly;g->Lz=Lz;g->dx=Lx/(double)nx;g->dy=Ly/(double)ny;g->dz=Lz/(double)nz;
    g->x=g->y=g->z=g->xc=g->yc=g->zc=g->cell_vol=NULL;return g;
}
void cfd_grid3d_destroy(CfdGrid3D*g){if(g){free(g->x);free(g->y);free(g->z);free(g->xc);free(g->yc);free(g->zc);free(g->cell_vol);free(g);}}

/* ===== Scalar Fields ===== */
CfdScalarField1D *cfd_scalar1d_create(int ncells,CfdFieldType ftype){
    if(ncells<1)return NULL;
    CfdScalarField1D*f=malloc(sizeof(CfdScalarField1D));if(!f)return NULL;
    f->ncells=ncells;f->ftype=ftype;f->loc=CFD_LOC_CENTER;f->name[0]=0;
    f->data=calloc((size_t)ncells,sizeof(double));if(!f->data){free(f);return NULL;}
    return f;
}
void cfd_scalar1d_destroy(CfdScalarField1D*f){if(f){free(f->data);free(f);}}
int cfd_scalar1d_fill(CfdScalarField1D*f,double v){if(!f||!f->data)return -1;for(int i=0;i<f->ncells;i++)f->data[i]=v;return 0;}

CfdScalarField2D *cfd_scalar2d_create(int nx,int ny,CfdFieldType ftype){
    if(nx<2||ny<2)return NULL;
    CfdScalarField2D*f=malloc(sizeof(CfdScalarField2D));if(!f)return NULL;
    f->nx=nx;f->ny=ny;f->ncells=nx*ny;f->ftype=ftype;f->loc=CFD_LOC_CENTER;f->name[0]=0;
    f->data=calloc((size_t)f->ncells,sizeof(double));if(!f->data){free(f);return NULL;}
    return f;
}
void cfd_scalar2d_destroy(CfdScalarField2D*f){if(f){free(f->data);free(f);}}
int cfd_scalar2d_copy(const CfdScalarField2D*src,CfdScalarField2D*dst){
    if(!src||!dst||!src->data||!dst->data)return -1;
    if(src->ncells!=dst->ncells)return -1;
    memcpy(dst->data,src->data,(size_t)src->ncells*sizeof(double));return 0;
}

/* ===== Vector Field ===== */
CfdVectorField2D *cfd_vector2d_create(int nx,int ny){
    if(nx<2||ny<2)return NULL;
    CfdVectorField2D*f=malloc(sizeof(CfdVectorField2D));if(!f)return NULL;
    f->nx=nx;f->ny=ny;f->ncells=nx*ny;f->u_loc=f->v_loc=f->p_loc=CFD_LOC_CENTER;f->name[0]=0;
    f->u=calloc((size_t)f->ncells,sizeof(double));
    f->v=calloc((size_t)f->ncells,sizeof(double));
    f->p=calloc((size_t)f->ncells,sizeof(double));
    if(!f->u||!f->v||!f->p){free(f->u);free(f->v);free(f->p);free(f);return NULL;}
    return f;
}
void cfd_vector2d_destroy(CfdVectorField2D*f){if(f){free(f->u);free(f->v);free(f->p);free(f);}}

/* ===== Coefficient Matrix ===== */
CfdMatrix2D *cfd_matrix2d_create(int nx,int ny){
    if(nx<1||ny<1)return NULL;
    CfdMatrix2D*m=malloc(sizeof(CfdMatrix2D));if(!m)return NULL;
    m->nx=nx;m->ny=ny;m->ncells=nx*ny;int n=m->ncells;
    m->aW=calloc((size_t)n,sizeof(double));m->aS=calloc((size_t)n,sizeof(double));
    m->aP=calloc((size_t)n,sizeof(double));m->aE=calloc((size_t)n,sizeof(double));
    m->aN=calloc((size_t)n,sizeof(double));m->b=calloc((size_t)n,sizeof(double));
    if(!m->aW||!m->aS||!m->aP||!m->aE||!m->aN||!m->b){cfd_matrix2d_destroy(m);return NULL;}
    return m;
}
void cfd_matrix2d_destroy(CfdMatrix2D*m){if(m){free(m->aW);free(m->aS);free(m->aP);free(m->aE);free(m->aN);free(m->b);free(m);}}
void cfd_matrix2d_zero(CfdMatrix2D*m){if(!m)return;int n=m->ncells;memset(m->aW,0,n*sizeof(double));memset(m->aS,0,n*sizeof(double));memset(m->aP,0,n*sizeof(double));memset(m->aE,0,n*sizeof(double));memset(m->aN,0,n*sizeof(double));memset(m->b,0,n*sizeof(double));}

/* ====================================================================
 * L2: Solver State
 * ==================================================================== */

void cfd_solver_state_init(CfdSolverState *state, CfdSolverType stype,
                            double tol, int max_iter)
{
    if (!state) return;
    state->solver_type = stype;
    state->tolerance = tol;
    state->max_iter = max_iter;
    state->n_iter = 0;
    state->initial_residual = 1.0;
    state->final_residual = 1.0;
    state->final_delta = 0.0;
    state->under_relaxation = 1.0;
    state->diverged = 0;
    state->history = NULL;
    state->history_capacity = 0;
    state->history_count = 0;
}

int cfd_solver_converged(const CfdSolverState *state)
{
    if (!state) return 0;
    return (state->final_residual < state->tolerance) && !state->diverged;
}

/* ====================================================================
 * L1: Boundary Condition Helpers
 * ==================================================================== */

int cfd_bc1d_default(CfdBCSet1D *bc)
{
    if (!bc) return -1;
    bc->left.type = CFD_BC_INLET_VELOCITY;
    bc->left.value = 0.0;
    bc->left.gradient = 0.0;
    bc->left.htc = 0.0;
    bc->left.T_ref = 300.0;
    bc->left.wall_velocity = 0.0;

    bc->right.type = CFD_BC_OUTLET_OUTFLOW;
    bc->right.value = 0.0;
    bc->right.gradient = 0.0;
    bc->right.htc = 0.0;
    bc->right.T_ref = 300.0;
    bc->right.wall_velocity = 0.0;
    return 0;
}

int cfd_bc2d_default(CfdBCSet2D *bc)
{
    if (!bc) return -1;
    /* West wall (no-slip) */
    bc->west.type = CFD_BC_WALL_NO_SLIP;
    bc->west.value = 0.0;
    bc->west.gradient = 0.0;
    bc->west.htc = 0.0;
    bc->west.T_ref = 300.0;
    bc->west.wall_velocity = 0.0;
    /* East outlet */
    bc->east.type = CFD_BC_OUTLET_OUTFLOW;
    bc->east.value = 0.0;
    bc->east.gradient = 0.0;
    bc->east.htc = 0.0;
    bc->east.T_ref = 300.0;
    bc->east.wall_velocity = 0.0;
    /* South wall */
    bc->south.type = CFD_BC_WALL_NO_SLIP;
    bc->south.value = 0.0;
    bc->south.gradient = 0.0;
    bc->south.htc = 0.0;
    bc->south.T_ref = 300.0;
    bc->south.wall_velocity = 0.0;
    /* North wall */
    bc->north.type = CFD_BC_WALL_NO_SLIP;
    bc->north.value = 0.0;
    bc->north.gradient = 0.0;
    bc->north.htc = 0.0;
    bc->north.T_ref = 300.0;
    bc->north.wall_velocity = 0.0;
    return 0;
}

/* ====================================================================
 * L3: Dimensionless Group Computation
 * ==================================================================== */

void cfd_fill_dimensionless_full(const CfdFluidProperty *fp,
                                  double U_ref, double L_ref,
                                  double dx_cell, double dt_step,
                                  CfdDimensionless *dless)
{
    if (!fp || !dless) return;

    double rho = fp->density;
    double mu  = fp->viscosity;
    double nu  = fp->kinematic_visc;
    double k   = fp->thermal_cond;
    double cp  = fp->specific_heat_cp;
    double alpha = k / (rho * cp);

    dless->Reynolds = rho * U_ref * L_ref / mu;
    dless->cell_Reynolds = rho * U_ref * dx_cell / mu;
    dless->Peclet = dless->Reynolds * fp->prandtl;
    dless->cell_Peclet = rho * U_ref * dx_cell / (k / cp);
    dless->Prandtl = fp->prandtl;
    dless->Schmidt = nu / CFD_SAFE_DIV(fp->mass_diffusivity, 1e-20);
    dless->Nusselt = 0.0; /* computed separately */
    dless->Sherwood = 0.0;
    dless->Grashof = 0.0; /* needs delta_T, g */;
    dless->Rayleigh = 0.0;
    dless->Courant = U_ref * dt_step / dx_cell;
    dless->Fourier_diff = alpha * dt_step / (dx_cell * dx_cell);
    dless->Stanton = CFD_SAFE_DIV(dless->Nusselt,
                                   dless->Reynolds * dless->Prandtl);
    dless->Mach = 0.0; /* needs speed_of_sound */
    dless->Weber = 0.0; /* needs surface_tension */
    dless->Froude = 0.0; /* needs g */
}

/* ====================================================================
 * L3: Fluid Property Initialization (Engineering Data)
 * ==================================================================== */

int cfd_fluid_air(CfdFluidProperty *fp)
{
    if (!fp) return -1;
    fp->density = CFD_AIR_DENSITY;
    fp->viscosity = CFD_AIR_VISCOSITY;
    fp->kinematic_visc = CFD_AIR_VISCOSITY / CFD_AIR_DENSITY;
    fp->thermal_cond = CFD_AIR_CONDUCTIVITY;
    fp->specific_heat_cp = 1007.0; /* J/(kg K) */
    fp->prandtl = CFD_AIR_PRANDTL;
    fp->thermal_expansion = 1.0 / 300.0; /* ideal gas, 1/T */
    fp->mass_diffusivity = 2.5e-5; /* m^2/s, typical gas */
    return 0;
}

int cfd_fluid_water(CfdFluidProperty *fp)
{
    if (!fp) return -1;
    fp->density = CFD_WATER_DENSITY;
    fp->viscosity = CFD_WATER_VISCOSITY;
    fp->kinematic_visc = CFD_WATER_VISCOSITY / CFD_WATER_DENSITY;
    fp->thermal_cond = CFD_WATER_CONDUCTIVITY;
    fp->specific_heat_cp = 4182.0;
    fp->prandtl = CFD_WATER_PRANDTL;
    fp->thermal_expansion = 2.07e-4; /* 1/K at 20C */
    fp->mass_diffusivity = 1.0e-9; /* m^2/s, typical liquid */
    return 0;
}

int cfd_fluid_oil(CfdFluidProperty *fp)
{
    if (!fp) return -1;
    fp->density = CFD_OIL_DENSITY;
    fp->viscosity = CFD_OIL_VISCOSITY;
    fp->kinematic_visc = CFD_OIL_VISCOSITY / CFD_OIL_DENSITY;
    fp->thermal_cond = CFD_OIL_CONDUCTIVITY;
    fp->specific_heat_cp = 1900.0;
    fp->prandtl = CFD_OIL_PRANDTL;
    fp->thermal_expansion = 7.0e-4;
    fp->mass_diffusivity = 1.0e-10;
    return 0;
}

int cfd_fluid_mercury(CfdFluidProperty *fp)
{
    if (!fp) return -1;
    fp->density = CFD_HG_DENSITY;
    fp->viscosity = CFD_HG_VISCOSITY;
    fp->kinematic_visc = CFD_HG_VISCOSITY / CFD_HG_DENSITY;
    fp->thermal_cond = CFD_HG_CONDUCTIVITY;
    fp->specific_heat_cp = 140.0;
    fp->prandtl = CFD_HG_PRANDTL;
    fp->thermal_expansion = 1.81e-4;
    fp->mass_diffusivity = 1.0e-9;
    return 0;
}

/* ====================================================================
 * L4: Residual Norms Computation
 * ==================================================================== */

void cfd_residual_norms_init(CfdResidualNorms *r)
{
    if (!r) return;
    r->L1_norm = 0.0;
    r->L2_norm = 0.0;
    r->L_inf_norm = 0.0;
    r->relative_L2 = 0.0;
    r->mass_imbalance = 0.0;
}

void cfd_compute_residual_norms(const double *residual, int n,
                                 const double *ref_residual,
                                 CfdResidualNorms *norms)
{
    if (!residual || n < 1 || !norms) return;
    double l1 = 0.0, l2 = 0.0, linf = 0.0;
    for (int i = 0; i < n; i++) {
        double r = fabs(residual[i]);
        l1 += r;
        l2 += r * r;
        if (r > linf) linf = r;
    }
    norms->L1_norm = l1;
    norms->L2_norm = sqrt(l2);
    norms->L_inf_norm = linf;

    if (ref_residual) {
        double ref_l2 = 0.0;
        for (int i = 0; i < n; i++) ref_l2 += ref_residual[i] * ref_residual[i];
        norms->relative_L2 = (ref_l2 > 0.0) ? sqrt(l2) / sqrt(ref_l2) : l2;
    } else {
        norms->relative_L2 = 0.0;
    }
    norms->mass_imbalance = 0.0;
}

/* ===== Additional Utilities ===== */
int cfd_grid1d_refine(const CfdGrid1D *coarse, CfdGrid1D **fine) {
    if (!coarse || !fine) return -1;
    *fine = cfd_grid1d_create(2 * coarse->nx, coarse->L);
    return (*fine) ? 0 : -1;
}
int cfd_grid2d_refine(const CfdGrid2D *coarse, CfdGrid2D **fine) {
    if (!coarse || !fine) return -1;
    *fine = cfd_grid2d_create(2*coarse->nx, 2*coarse->ny, coarse->Lx, coarse->Ly);
    if (!(*fine)) return -1;
    return cfd_grid2d_allocate_all(*fine);
}
double cfd_grid1d_min_spacing(const CfdGrid1D *grid) { return grid ? grid->dx : 0.0; }
double cfd_grid2d_min_spacing(const CfdGrid2D *grid) {
    return (grid && grid->dx > 0 && grid->dy > 0) ?
        CFD_MIN(grid->dx, grid->dy) : 0.0;
}
int cfd_grid2d_cell_count(const CfdGrid2D *grid) { return grid ? grid->ncells : 0; }
int cfd_grid1d_cell_count(const CfdGrid1D *grid) { return grid ? grid->nx : 0; }

void cfd_vector2d_set_uniform(CfdVectorField2D *vel, double u0, double v0, double p0) {
    if (!vel) return;
    for (int i = 0; i < vel->ncells; i++) { vel->u[i] = u0; vel->v[i] = v0; vel->p[i] = p0; }
}
double cfd_vector2d_max_velocity(const CfdVectorField2D *vel) {
    if (!vel) return 0.0;
    double max_v = 0.0;
    for (int i = 0; i < vel->ncells; i++) {
        double vmag = sqrt(vel->u[i]*vel->u[i] + vel->v[i]*vel->v[i]);
        if (vmag > max_v) max_v = vmag;
    }
    return max_v;
}

double cfd_field_min(const CfdScalarField1D *field) {
    if (!field || field->ncells < 1) return 0.0;
    double v = field->data[0];
    for (int i = 1; i < field->ncells; i++) if (field->data[i] < v) v = field->data[i];
    return v;
}
double cfd_field_max(const CfdScalarField1D *field) {
    if (!field || field->ncells < 1) return 0.0;
    double v = field->data[0];
    for (int i = 1; i < field->ncells; i++) if (field->data[i] > v) v = field->data[i];
    return v;
}
double cfd_field2d_min(const CfdScalarField2D *field) {
    if (!field || field->ncells < 1) return 0.0;
    double v = field->data[0];
    for (int i = 1; i < field->ncells; i++) if (field->data[i] < v) v = field->data[i];
    return v;
}
double cfd_field2d_max(const CfdScalarField2D *field) {
    if (!field || field->ncells < 1) return 0.0;
    double v = field->data[0];
    for (int i = 1; i < field->ncells; i++) if (field->data[i] > v) v = field->data[i];
    return v;
}

double cfd_field_average(const CfdScalarField1D *field) {
    if (!field || field->ncells < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < field->ncells; i++) sum += field->data[i];
    return sum / (double)field->ncells;
}
double cfd_field_rms(const CfdScalarField1D *field) {
    if (!field || field->ncells < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < field->ncells; i++) sum += field->data[i] * field->data[i];
    return sqrt(sum / (double)field->ncells);
}

void cfd_solver_state_record(CfdSolverState *state, double residual) {
    if (!state) return;
    if (state->history_count < state->history_capacity && state->history) {
        state->history[state->history_count].iter = state->history_count;
        state->history[state->history_count].residual = residual;
        state->history_count++;
    }
}
int cfd_convergence_criterion_check(const CfdConvergenceCriterion *c,
                                      double residual, double initial_residual, int iter) {
    if (!c) return 0;
    if (iter < c->min_iter) return 0;
    if (residual > c->div_tolerance) return 0;
    if (residual < c->abs_tolerance) return 1;
    if (initial_residual > 0 && residual/initial_residual < c->rel_tolerance) return 1;
    return 0;
}
