/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing authors:
     L.T. Kong, C. Campana, C. Denniston and M. Muser

   Contact:
     Department of Applied Mathematics, University of Western Ontario
     London, ON, Canada N6A 5B7

     mmuser@uwo.ca, cdennist@uwo.ca, konglt@gmail.com
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "fix_wall_rough.h"
#include "fft3d_wrap.h"
#include "random_park.h"
#include "memory.h"
#include "comm.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 * syntax:
 * fix ID g-ID wall/rough style nx ny seed qs ql hurst rmss pottype args keyword/values
 *     stype = xlo or xhi or ylo or yhi or zlo or zhi
 *     nx,ny = discretization in the plane direction; must >1
 *     seed  = random number seed to generate wall
 *     qs,ql = low/high limit of wave, in unit of 2pi/a, [0,1]
 *     hurst = Hurst exponent of rough wall
 *     rmss  = desired RMS slope of wall
 *     pottype = potential of wall-particle interaction
 *          hardwall = None; simple reset particle at wall position if pass through
 *          morse    = d0 alpha r0 cutoff; Morse potential
 *          lj       = epsilon sigma cutoff; LJ 12-6
 *          lj93     = epsilon sigma cutoff; LJ 9-3
 *          repulsive= d0 sigma cutoff; Exponentially repulsive potential
 *     keyword = wmax or wmin or wmean or file
 *       wmax = max_value; shift topmost position of wall to wmax
 *       wmin = min_value; shift lowest position of wall to wmin
 *       wmean = mean_value; shift mean position to wmean
 *       file  = file name; output wall info to file
 * ---------------------------------------------------------------------*/

FixWallRough::FixWallRough(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (domain->triclinic != 0) error->all(FLERR,"Fix wall/rough works with orthogonal lattice only");
  if (narg < 12) error->all(FLERR,"Illegal fix wall/rough command");

  // wall position
  if (strcmp(arg[3],"xlo") == 0) {
    zdim = 0; xdim=1; ydim=2;
    side = -1;
  } else if (strcmp(arg[3],"xhi") == 0) {
    zdim = 0; xdim=1; ydim=2;
    side = 1;
  } else if (strcmp(arg[3],"ylo") == 0) {
    zdim = 1; xdim=0; ydim=2;
    side = -1;
  } else if (strcmp(arg[3],"yhi") == 0) {
    zdim = 1; xdim=0; ydim=2;
    side = 1;
  } else if (strcmp(arg[3],"zlo") == 0) {
    zdim = 2; xdim=0; ydim=1;
    side = -1;
  } else if (strcmp(arg[3],"zhi") == 0) {
    zdim = 2; xdim=0; ydim=1;
    side = 1;
  } else error->all(FLERR,"Illegal fix wall/rough command");

  if (domain->periodicity[zdim]) error->all(FLERR,"Cannot use wall in periodic dimension");

  // wall generation controling parameters
  int iarg=4;
  nx   = atoi(arg[iarg++]);
  ny   = atoi(arg[iarg++]);
  seed = atoi(arg[iarg++]);
  qs   = atof(arg[iarg++]);
  ql   = atof(arg[iarg++]);
  hurst= atof(arg[iarg++]);
  rmss = atof(arg[iarg++]);
  
  if (nx<1||ny<1||seed<1||qs<0.||qs>1.||ql<0.||ql>1.||ql<qs
     ||hurst<0.||rmss<0.) error->all(FLERR,"Illegal fix wall/rough command");
  
  // wall-particle interaction type
  if (strcmp(arg[iarg],"hardwall") == 0) {
    pottype = 0;

  } else if (strcmp(arg[iarg],"morse") == 0) {
    if (iarg+5>narg) error->all(FLERR,"Insufficient command line option for fix wall/rough!");
    pottype = 1;
    d0      = atof(arg[++iarg]);
    alpha   = atof(arg[++iarg]);
    r0      = atof(arg[++iarg]);
    cutoff  = atof(arg[++iarg]);
    
  } else if (strcmp(arg[iarg],"lj") == 0) {
    if (iarg+4>narg) error->all(FLERR,"Insufficient command line option for fix wall/rough!");
    pottype = 2;
    epsilon = atof(arg[++iarg]);
    sigma   = atof(arg[++iarg]);
    cutoff  = atof(arg[++iarg]);

  } else if (strcmp(arg[iarg],"lj93") == 0){
    if (iarg+4>narg) error->all(FLERR,"Insufficient command line option for fix wall/rough!");
    pottype = 3;
    epsilon = atof(arg[++iarg]);
    sigma   = atof(arg[++iarg]);
    cutoff  = atof(arg[++iarg]);

  } else if (strcmp(arg[iarg],"repulsive") ==0){
    if (iarg+4>narg) error->all(FLERR,"Insufficient command line option for fix wall/rough!");
    pottype = 4;
    d0      = atof(arg[++iarg]);
    sigma   = atof(arg[++iarg]);
    cutoff  = atof(arg[++iarg]);

    rsigma  = 1./sigma;

  } else error->all(FLERR,"Illegal fix wall/rough command");
  iarg++;
 
  int shift; // default position of rough wall
  if (side == -1){ // on bottom, the maximum sits at the low limit of simulation box
    shift = 1;
    wallshift = domain->boxlo[zdim];
  } else {         // on top, the minimum sits at the high limit of simulation box
    shift = 2;
    wallshift = domain->boxhi[zdim];
  }
  fname  = NULL;

  while (iarg+1<narg){ // read keyword
    if (strcmp(arg[iarg],"wmax") == 0){ // top limit of wall
      if (iarg+2>narg) error->all(FLERR,"Insufficient command line option for fix wall/rough!");
      shift = 1;
      wallshift = atof(arg[++iarg]);

    } else if (strcmp(arg[iarg],"wmin") == 0){ // bottom limit of wall
      if (iarg+2>narg) error->all(FLERR,"Insufficient command line option for fix wall/rough!");
      shift = 2;
      wallshift = atof(arg[++iarg]);

    } else if (strcmp(arg[iarg],"wmean") == 0){ // mean position of wall
      if (iarg+2>narg) error->all(FLERR,"Insufficient command line option for fix wall/rough!");
      shift = 3;
      wallshift = atof(arg[++iarg]);

    } else if (strcmp(arg[iarg],"file") ==0){ // output wall info
      if (iarg+2>narg) error->all(FLERR,"Insufficient command line option for fix wall/rough!");
      if (fname) delete []fname;
      int n = strlen(arg[++iarg])+1;
      fname = new char[n];
      strcpy(fname,arg[iarg]);

    } else error->all(FLERR,"Illegal fix wall/rough command");

    iarg++;
  }

  if (pottype > 0){ // hard-wall does not compute energy/force; the others do
    scalar_flag = 1;
    vector_flag = 1;
    size_vector = 3;
    global_freq = 1;
    extscalar = 1;
    extvector = 1;
  }

  wall_flag = 0;
  wall[0] = wall[1] = wall[2] = wall[3] = 0.0; // energy and f on wall

  CreateSurface(shift);
}

/* ---------------------------------------------------------------------- */

FixWallRough::~FixWallRough()
{
  memory->destroy(h);
  memory->destroy(hx);
  memory->destroy(hy);
  memory->destroy(hxy);

  if (fname) delete []fname;
}

/* ---------------------------------------------------------------------- */

int FixWallRough::setmask()
{
  int mask = 0;
  if (pottype == 0)
    mask |= FINAL_INTEGRATE;
  else {
    mask |= POST_FORCE;
    mask |= THERMO_ENERGY;
    mask |= POST_FORCE_RESPA;
    mask |= MIN_POST_FORCE;
  }
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallRough::init()
{
  // initialize potential offset and other coeff's
  if (pottype == 1){
    ald0 = 2.*alpha*d0;
    double armr0 = -alpha*(cutoff-r0);
    double exp_tmp = exp(armr0);
    offset = d0*exp_tmp*(exp_tmp-2.);

  } else if (pottype == 2){
    coeff1 = 48. * epsilon * pow(sigma,12.0);
    coeff2 = 24. * epsilon * pow(sigma,6.0);
    coeff3 = 4.0 * epsilon * pow(sigma,12.0);
    coeff4 = 4.0 * epsilon * pow(sigma,6.0);

    double r2inv = 1.0/(cutoff*cutoff);
    double r6inv = r2inv*r2inv*r2inv;
    offset = r6inv*(coeff3*r6inv - coeff4);

  } else if (pottype == 3){
    coeff1 = 6.0/5.0 * epsilon * pow(sigma,9.0);
    coeff2 = 3.0 * epsilon * pow(sigma,3.0);
    coeff3 = 2.0/15.0 * epsilon * pow(sigma,9.0);
    coeff4 = epsilon * pow(sigma,3.0);
  
    double rinv = 1.0/cutoff;
    double r2inv = rinv*rinv;
    double r4inv = r2inv*r2inv;
    offset = coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv;
  } else if (pottype == 4){
    offset = d0*exp(-cutoff*rsigma);
  }

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallRough::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRough::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */
double FixWallRough::memory_usage()
{
  double bytes = sizeof(double)*nx*ny*4;

  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixWallRough::post_force(int vflag)
{
  wall_flag = 0;
  wall[0] = wall[1] = wall[2] = wall[3] = 0.0;

  if      (pottype == 1){Morse();}
  else if (pottype == 2){LJ126();}
  else if (pottype == 3){LJ93();}
  else if (pottype == 4){Repulsive();}
}

/* ---------------------------------------------------------------------
 * private method, to calculate the force and energy of LJ 12-6 
 * interaction between rough wall and particle.
 * --------------------------------------------------------------------*/

void FixWallRough::LJ126()
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask  = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double delta,rinv,r2inv,r6inv,fwall;
  double xpos,ypos,zpos,wpos;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      domain->remap(x[i], image[i]);
      xpos = x[i][xdim]-domain->boxlo[xdim];
      ypos = x[i][ydim]-domain->boxlo[ydim];
      zpos = x[i][zdim];
      wpos = bicuint(xpos,ypos);

      delta = fabs(wpos-zpos);
      // because of rough wall, particle might goes under wall,
      // so always make the force on particle to away from wall
      if (delta>=cutoff) continue;
      
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r6inv = r2inv*r2inv*r2inv;
      fwall = double(side) * r6inv*(coeff1*r6inv - coeff2) * rinv;
      f[i][zdim]   -= fwall;
      wall[0]      += r6inv*(coeff3*r6inv - coeff4) - offset;
      wall[zdim+1] += fwall;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRough::LJ93()
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask  = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double delta,rinv,r2inv,r4inv,r10inv,fwall;
  double xpos, ypos, zpos, wpos;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      domain->remap(x[i], image[i]);
      xpos = x[i][xdim]-domain->boxlo[xdim];
      ypos = x[i][ydim]-domain->boxlo[ydim];
      zpos = x[i][zdim];
      wpos = bicuint(xpos,ypos);

      delta = fabs(wpos-zpos);
      if (delta>=cutoff) continue;

      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r4inv = r2inv*r2inv;
      r10inv = r4inv*r4inv*r2inv;
      fwall = (coeff1*r10inv - coeff2*r4inv) * double(side);
      f[i][zdim]   -= fwall;
      wall[0]      += coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv - offset;
      wall[zdim+1] += fwall;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRough::Morse()
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask  = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double xpos,ypos,zpos,wpos;
  double delta, exp_tmp, fwall;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      domain->remap(x[i], image[i]);
      xpos = x[i][xdim]-domain->boxlo[xdim];
      ypos = x[i][ydim]-domain->boxlo[ydim];
      zpos = x[i][zdim];
      wpos = bicuint(xpos,ypos);

      delta = fabs(wpos-zpos);
      if (delta>=cutoff) continue;

      exp_tmp = exp(-alpha*(delta-r0));
      fwall = double(side)*ald0*exp_tmp*(exp_tmp-1.);
      
      f[i][zdim]   -= fwall;
      wall[0]      += d0*exp_tmp*(exp_tmp-2.)-offset;
      wall[zdim+1] += fwall;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRough::Repulsive()
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask  = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double xpos,ypos,zpos,wpos;
  double delta, dexp, fwall;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      domain->remap(x[i], image[i]);
      xpos = x[i][xdim]-domain->boxlo[xdim];
      ypos = x[i][ydim]-domain->boxlo[ydim];
      zpos = x[i][zdim];
      wpos = bicuint(xpos,ypos);

      delta = (wpos-zpos)*side;
      if (delta>=cutoff) continue;

      dexp = d0*exp(-delta*rsigma);
      fwall = rsigma*dexp*side;

      f[i][zdim]   -= fwall;
      wall[0]      += dexp - offset;
      wall[zdim+1] += fwall;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRough::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRough::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
 * hard-wall, simply reset the position of atom as the wall posistion and
 * zero its corresponding velocity if that atom goes through the wall.
 * --------------------------------------------------------------------*/

void FixWallRough::final_integrate()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask  = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double xpos, ypos, zpos, wpos;

  for (int i=0; i<nlocal; i++){
    if (mask[i] & groupbit){
      domain->remap(x[i], image[i]);
      xpos = x[i][xdim]-domain->boxlo[xdim];
      ypos = x[i][ydim]-domain->boxlo[ydim];
      zpos = x[i][zdim];
      wpos = bicuint(xpos, ypos);
      if (double(side)*(wpos-zpos)<0.){x[i][zdim]=wpos; v[i][zdim]=0.;}
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRough::CreateSurface(int flag)
{
  double **g,**p, tpi = 8.*atan(1.);
  RanPark *random;

  dx  = domain->h[xdim]/nx;
  dy  = domain->h[ydim]/ny;
  rdx = 1./dx;
  rdy = 1./dy;
  qs *= MIN(rdx,rdy);
  ql *= MAX(rdx,rdy);

  g   = memory->create(g, nx,ny,"fix_wall/rough:g");
  p   = memory->create(p, nx,ny,"fix_wall/rough:p");
  h   = memory->create(h, nx,ny,"fix_wall/rough:h");
  
  random = new RanPark(lmp,seed);

  double qx, qy, qr;
  double expon = -hurst-(domain->dimension-1.)*0.5;
  for (int i=0; i<nx; i++){
    qx = double(i)/nx*rdx;
    for (int j=0; j<ny; j++){
      qy = double(j)/ny*rdy;
      qr = sqrt(qx*qx+qy*qy);
      
      if (qr>=qs && qr <=ql){
        g[i][j] = random->gaussian() * pow(qr,expon);
        p[i][j] = random->uniform()  * tpi;
      } else {
        g[i][j] = p[i][j] = 0.;
      }
    }
  }
  g[0][0] = 0.; 
  GetHviaFFT(g,p); 

  delete random;
  memory->destroy(g);
  memory->destroy(p);

  // rescale height to get desired RMS slope
  double scale = rmss/RMSSlope();
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++) h[i][j] *= scale;
  }

  // re-positioning the wall to desired z position
  double wmax=h[0][0], wmin=h[0][0], wmean=0., shift=0.;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      wmax   = MAX(wmax,h[i][j]);
      wmin   = MIN(wmin,h[i][j]);
      wmean += h[i][j];
    }
  }
  wmean /= double(nx*ny);
  if      (flag == 1) shift = wmax;
  else if (flag == 2) shift = wmin;
  else if (flag == 3) shift = wmean;
  shift -= wallshift;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++) h[i][j] -= shift;
  }
  wmax  -= shift;
  wmin  -= shift;
  wmean -= shift;

  // output wall info if required
  if (fname && comm->me==0){
    FILE *fp;
    char dir[5];
    dir[0] ='x';dir[1]='y';dir[2]='z';dir[3]='L';dir[4]='H';
    fp = fopen(fname,"w");
    fprintf(fp,"# Rough wall of RMS slope= %lg with dimension %d x %d\n",rmss,nx,ny);
    fprintf(fp,"# at %c%c ", dir[zdim],dir[3+(side+1)/2]);
    fprintf(fp,", dx= %lg , dy= %lg\n",dx,dy);
    fprintf(fp,"# Wmax= %lg , Wmin= %lg , Wmean= %lg\n",wmax,wmin,wmean);
    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
        fprintf(fp,"%lg %lg %lg\n",double(i)*dx,double(j)*dy,h[i][j]);
      }
    }
    fclose(fp);
  }

  // compute differencials for bicubic interpolation
  hx  = memory->create(hx,  nx,ny,"fix_wall/rough:hx");
  hy  = memory->create(hy,  nx,ny,"fix_wall/rough:hy");
  hxy = memory->create(hxy, nx,ny,"fix_wall/rough:hxy");

  double r2dx=rdx*0.5, r2dy=rdy*0.5, r2dxdy = r2dx*r2dy;
  for (int i=0; i<nx; i++){
    int xp1=(i+1)%nx;
    int xm1=(i-1+nx)%nx;
    for (int j=0; j<ny; j++){
      int yp1=(j+1)%ny;
      int ym1=(j-1+ny)%ny;
      hx[i][j]  = (h[xp1][j]-h[xm1][j])*r2dx;
      hy[i][j]  = (h[i][yp1]-h[j][ym1])*r2dy;
      hxy[i][j] = (h[xp1][yp1]-h[xp1][ym1]-h[xm1][yp1]+h[xm1][ym1])*r2dxdy;
    }
  }

  return;
}

/* ----------------------------------------------------------------------
 * private method to calculate the root-mean-square slope of wall,
 * periodicity is made use of.
 * --------------------------------------------------------------------*/

double FixWallRough::RMSSlope()
{
  double sumx=0., sumy=0., delx, dely;
  for (int i=0; i<nx; i++){
    int im1 = (i-1+nx)%nx;
    for (int j=0; j<ny; j++){
      int jm1 = (j-1+ny)%ny;
      delx = h[i][j] - h[im1][j];
      dely = h[i][j] - h[i][jm1];
      sumx += delx*delx;
      sumy += dely*dely;
    }
  }
  return sqrt((sumx*(rdx*rdx)+sumy*(rdy*rdy))/double(nx*ny));
}

/* ----------------------------------------------------------------------
 * energy of wall interaction
 * --------------------------------------------------------------------*/

double FixWallRough::compute_scalar()
{
  // only sum across procs one time
  if (wall_flag == 0) {
    MPI_Allreduce(wall,wall_all,4,MPI_DOUBLE,MPI_SUM,world);
    wall_flag = 1;
  }
  return wall_all[0];
}

/* ----------------------------------------------------------------------
 * components of force on wall
 * --------------------------------------------------------------------*/

double FixWallRough::compute_vector(int n)
{
  // only sum across procs one time
  if (wall_flag == 0) {
    MPI_Allreduce(wall,wall_all,4,MPI_DOUBLE,MPI_SUM,world);
    wall_flag = 1;
  }
  return wall_all[n+1];
}

/* ----------------------------------------------------------------------
 * private method, to get the coefficients for bicubic interpolation
 * --------------------------------------------------------------------*/

void FixWallRough::bcucof(double *y,  double *y1, double *y2,
                          double *y12,double  d1, double  d2, double *c)
{
  int i, j;
  double d1d2, x[16],xx;
  const double wt[][16] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1};
  /*-------------------------------------------------------------------*/
  d1d2 = d1 * d2;
  for (i=0; i<4; i++){
    x[i]   = y[i];
    x[i+4] = y1[i] * d1;
    x[i+8] = y2[i] * d2;
    x[i+12]= y12[i]* d1d2;
  }
  for (i=0; i<16; i++){
    xx = 0.;
    for (j=0; j<16; j++) xx += wt[i][j] * x[j];
    c[i] = xx;
  }
  return;
}

/* ----------------------------------------------------------------------
 * private method, to get the wall coordinate at the given plane position
 * by bi-cubic interpolation
 * --------------------------------------------------------------------*/

double FixWallRough::bicuint(double x, double y)
{
  double y0[4],y1[4],y2[4],y12[4], c[4][4];
  double t = x*rdx, u = y*rdy, res = 0.;
  int i = int(t)%nx, j = int(u)%ny;
  int xp1 = (i+1)%nx, yp1 = (j+1)%ny;

  y0[0] = h[i][j];
  y0[1] = h[xp1][j];
  y0[2] = h[xp1][yp1];
  y0[3] = h[i][yp1];

  y1[0] = hx[i][j];
  y1[1] = hx[xp1][j];
  y1[2] = hx[xp1][yp1];
  y1[3] = hx[i][yp1];

  y2[0] = hy[i][j];
  y2[1] = hy[xp1][j];
  y2[2] = hy[xp1][yp1];
  y2[3] = hy[i][yp1];

  y12[0] = hxy[i][j];
  y12[1] = hxy[xp1][j];
  y12[2] = hxy[xp1][yp1];
  y12[3] = hxy[i][yp1];

  bcucof(y0,y1,y2,y12,dx,dy,c[0]);

  t -= floor(t); u -= floor(u);
  for (i=3; i>=0; i--) res = t * res + (c[i][3]*u+c[i][2]*u+c[i][1])*u + c[i][0];

  return res;
}

/* ----------------------------------------------------------------------
 * private method to get the wall position via FFT
 * --------------------------------------------------------------------*/

void FixWallRough::GetHviaFFT(double **g, double **p)
{
  FFT3d *fft;
  int nxlo, nxhi, mynpt, mysize;
  int *fft_cnts,*fft_disp,*locnpt;
  int me = comm->me, nprocs = comm->nprocs;
  double *buf, *fftdata;

  locnpt = new int [nprocs];
  for (int iproc=0; iproc<nprocs; iproc++){
    locnpt[iproc] = nx/nprocs;
    if (iproc >= (nprocs-nx%nprocs)) locnpt[iproc] += 1;
  }
  nxlo = 0;
  for (int i=0; i<comm->me; i++) nxlo += locnpt[i];
  nxhi  = nxlo + locnpt[me]-1;
  mynpt = locnpt[me]*ny;
  fft_cnts = new int[nprocs];
  fft_disp = new int[nprocs];
  fft_disp[0] = 0;
  for (int i=0; i<nprocs; i++) fft_cnts[i] = locnpt[i]*ny;
  for (int i=1; i<nprocs; i++) fft_disp[i] = fft_disp[i-1] + fft_cnts[i-1];

  fftdata = (double *) memory->smalloc(MAX(1,mynpt)*2*sizeof(double),"FixWallRough:fftdata");
  buf     = (double *) memory->smalloc(MAX(1,mynpt)  *sizeof(double),"FixWallRough:buf");
  fft     = new FFT3d(lmp,world,1,ny,nx,0,0,0,ny-1,nxlo,nxhi,0,0,0,ny-1,nxlo,nxhi,0,0,&mysize);

  // to fill data and do FFT
  int m = 0;
  for (int i=nxlo; i<nxhi; i++){
    for (int j=0; j<ny; j++){
      fftdata[m++] = g[i][j]*cos(p[i][j]);
      fftdata[m++] = g[i][j]*sin(p[i][j]);
    }
  }
  fft->compute(fftdata, fftdata, 1);
  m = 0;
  for (int i=0; i<mynpt; i++){ buf[i] = fftdata[m]; m += 2; }

  // Gather data 
  MPI_Allgatherv(buf,mynpt,MPI_DOUBLE,h[0],fft_cnts,fft_disp,MPI_DOUBLE,world);

  // free allocated memories
  memory->sfree(buf);
  memory->sfree(fftdata);
  delete fft;
  delete []fft_cnts;
  delete []fft_disp;
  delete []locnpt;
  
  return;
}
/* -------------------------------------------------------------------- */
