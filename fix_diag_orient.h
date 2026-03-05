/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------
   Customized fix to diagnose the orientation of molecules based on a given
   bond, angle, or dihedral type.
   The vector that connects atoms for 1-2 of bond, 1-3 of angle, or 1-4 of
   dihedral will be used to do the measurement and output the tensor.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(diag/orient,FixDiagOrient);
// clang-format on
#else

#ifndef LMP_FIX_DIAG_ORIENT_H
#define LMP_FIX_DIAG_ORIENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDiagOrient : public Fix {
 public:
  FixDiagOrient(class LAMMPS *, int, char **);
  ~FixDiagOrient() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  double compute_vector(int) override;
  double memory_usage() override;  // not used before

 private:
  int me, nprocs, first;
  int id_type;
  int noutput;
  int checkchg;
  double qref;

  char *base, *prefix;

  int Nbx,Nby,Nbz,restartflg;

  double **altogether;

  int ncount, num_ave;
  int subNbx,subNby,subNbz;
  double sublosft[3],subsize[3];

  // to enable vector output
  int cv_flag;
  double vector_me[7], vector_all[7];

};

}    // namespace LAMMPS_NS

#endif
#endif
