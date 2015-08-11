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

#ifdef FIX_CLASS

FixStyle(wall/rough,FixWallRough)

#else

#ifndef FIX_WALL_ROUGH_H
#define FIX_WALL_ROUGH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallRough : public Fix {
 public:
  FixWallRough(class LAMMPS *, int, char **);
  ~FixWallRough();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  void final_integrate();
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  int side,xdim,ydim,zdim;
  int nx,ny,seed,pottype;
  double hurst,rmss,qs,ql;
  double dx,dy,rdx,rdy,wallshift;
  double **h,**hx,**hy,**hxy;
  char *fname;

  int wall_flag;
  double wall[4],wall_all[4];
  int nlevels_respa;

  double epsilon,sigma,cutoff; // For LJ interaction
  double coeff1,coeff2,coeff3,coeff4,offset;

  double d0,alpha,r0,ald0;     // For Morse interaction: d0, alpha, r0, alpha*d0
  double rsigma;               // For exponential potential: d0, sigma

  void Morse();
  void LJ126();
  void LJ93();
  void Repulsive();

  void CreateSurface(int);
  void GetHviaFFT(double**,double**);
  double RMSSlope();
  void bcucof(double*,double*,double*,double*,double,double,double*);
  double bicuint(double,double);
  
};
}
#endif
#endif
