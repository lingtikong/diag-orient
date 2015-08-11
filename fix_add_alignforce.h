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

FixStyle(alignforce,FixAddAlignForce)

#else


#ifndef FIX_ADD_ALIGNFORCE_H
#define FIX_ADD_ALIGNFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAddAlignForce : public Fix {
 public:
  FixAddAlignForce(class LAMMPS *, int, char **);
  ~FixAddAlignForce();
  int setmask();
  void init();
  void post_force(int);

 private:
  int    nadd;
  int    *fdir;
  int    *btype;
  double *f2add;
  double *qsel;
  double zlo1, zhi1, zlo2, zhi2;
};

}

#endif
#endif
