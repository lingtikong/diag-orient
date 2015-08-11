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
#include "fix_add_alignforce.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "math.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

/* To add force on atom with perticular charge in a perticular bond; if this
   atom lies to the left of the other atom in the bond, a negative force is
   added, otherwise a positive force is added. */

FixAddAlignForce::FixAddAlignForce(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix addalignforce command");

  scalar_flag = 0;
  vector_flag = 0;

  nadd   = (narg-3)/4;
  if (nadd < 1) error->all(FLERR,"No align force necessary at all!");
  
  fdir   = new int[nadd];
  btype  = new int[nadd];
  f2add  = new double[nadd];
  qsel   = new double[nadd];
  
  int n = 2;
  for (int i=0; i<nadd; i++){
    fdir[i]  = atoi(arg[n+1]); // direction to add align force
    f2add[i] = atof(arg[n+2]); // value of force to add, in force unit
    btype[i] = atoi(arg[n+3]); // bond type 1
    qsel[i]  = atof(arg[n+4]); // force will be added to atom has bond type 1 with this charge
    if ( fdir[i] < 0 || fdir[i] > 2 ) error->all(FLERR,"Wrong direction to add force!");
    n+=4;
  }

}

/* ---------------------------------------------------------------------- */

FixAddAlignForce::~FixAddAlignForce()
{
  delete []fdir;
  delete []btype;
  delete []f2add;
  delete []qsel;
}

/* ---------------------------------------------------------------------- */

void FixAddAlignForce::init()
{
  zlo1 = domain->boxlo[2] + 0.10*domain->zprd;
  zhi1 = domain->boxlo[2] + 0.45*domain->zprd;
  zlo2 = domain->boxlo[2] + 0.60*domain->zprd;
  zhi2 = domain->boxlo[2] + 0.95*domain->zprd;
}


/* ---------------------------------------------------------------------- */

int FixAddAlignForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddAlignForce::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  double *q  = atom->q;
  double const zero = 1.e-4;
  int **bondlist = neighbor->bondlist;
  int  nbondlist = neighbor->nbondlist;
  int i1, i2, zpos, left, isel, ffac;
  double z1, z2, zb;

  for (int n = 0; n < nbondlist; n++){
    for (int i=0; i < nadd; i++){
      if (bondlist[n][2] == btype[i]){
        i1 = bondlist[n][0];
        i2 = bondlist[n][1];
        z1 = x[i1][2];
        z2 = x[i2][2];

        while ( z1 > domain->boxhi[2] ) z1 -= domain->zprd;
        while ( z1 < domain->boxlo[2] ) z1 += domain->zprd;
        while ( z2 > domain->boxhi[2] ) z2 -= domain->zprd;
        while ( z2 < domain->boxlo[2] ) z2 += domain->zprd;

        zb = (z1+z2)*0.5;
        zpos = 0;
        if ( zb > zlo1 && zb < zhi1 ){zpos = 1;}
        else if ( zb > zlo2 && zb < zhi2){zpos = -1;}
        if (zpos == 0) break;

        left = -1;
        if (x[i1][0] > x[i2][0]) left = 1;

        isel = -1;
        if ( fabs(q[i1]-qsel[i]) < zero){isel = i1; ffac = left;}
        else if ( fabs(q[i2]-qsel[i]) < zero){isel = i2; ffac=-left;}
        if (isel < 0) continue;
        
        int cdir = fdir[i];
        if ( cdir == 2) ffac *= zpos;
        f[isel][cdir] += double(ffac)*f2add[i];
      }
    }
  }
}
