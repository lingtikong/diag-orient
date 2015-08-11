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

FixStyle(orient,DiagOrientation)

#else

#ifndef ORIENTATION_H
#define ORIENTATION_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class DiagOrientation : public Fix {
 public:
  DiagOrientation(class LAMMPS *, int, char **);
  ~DiagOrientation();
  int setmask();
  void init();
  void end_of_step();
  double compute_vector(int);

 private:
  int me, nprocs, first;
  int orient_base, id_type;
  int noutput;
  int checkchg;
  double qref;

  char *prefix;

  int Nbx,Nby,Nbz,restartflg;

  double (*altogether)[10];

  int ncount, num_ave;
  int subNbx,subNby,subNbz;
  double sublosft[3],subsize[3];
  // to enable vector output
  int cv_flag;
  double vector_me[7], vector_all[7];
};

}
#endif
#endif
