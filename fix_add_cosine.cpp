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

#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "fix_add_cosine.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixAddCosine::FixAddCosine(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix addcosine command");
  if (domain->nonperiodic==1) error->all(FLERR,"Fix addcosine requires orthogonal box");

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  numadd = narg/3-1;
  dims   = new int[numadd];
  deltas = new double[numadd];
  ks     = new double[numadd];
  

  int iarg=2, idum=0;
  while ( ++iarg < narg ){
    if (strcmp(arg[iarg],"x") == 0) dims[idum] = 0;
    else if (strcmp(arg[iarg],"y") == 0) dims[idum] = 1;
    else if (strcmp(arg[iarg],"z") == 0) dims[idum] = 2;
    else error->all(FLERR,"Illegal fix addcosine command");

    deltas[idum] = atof(arg[++iarg]);
    ks[idum]     = atof(arg[++iarg]);
    idum++;
  }
  if (idum != numadd) error->all(FLERR,"Illegal fix addcosine command");

  for (int i=0; i<numadd; i++){
    int idim = dims[i];
    ks[i] *= 8.*atan(1.)/domain->prd[idim];
  }
  
  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
}

/* ---------------------------------------------------------------------- */
FixAddCosine::~FixAddCosine()
{
  delete []dims;
  delete []deltas;
  delete []ks;
}
/* ---------------------------------------------------------------------- */

int FixAddCosine::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddCosine::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAddCosine::setup(int vflag)
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

void FixAddCosine::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddCosine::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      foriginal[0] += f[i][0];
      foriginal[1] += f[i][1];
      foriginal[2] += f[i][2];
      for (int j=0; j<numadd; j++){
        int idim = dims[j];
        f[i][idim] += deltas[j] * cos(ks[j]*atom->x[i][idim]);
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixAddCosine::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddCosine::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAddCosine::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n];
}
