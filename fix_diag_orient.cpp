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

#include "fix_diag_orient.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "memory.h"
#include "mpi.h"

#include "string.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

DiagOrientation::DiagOrientation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // sample inputfile call: //
  // fix  ID g-ID orient nevery orient_base type_id prefix noutput xbin ybin zbin flag [qsum]
  // where
  // ID, g-ID is as other fix
  // nevery is the frequency to measure orientation
  // orient_base indicate the orientation is measured according to a bond (2),a angle (3) or a dihedral (4)
  // type_id is the id of the bond type or angle type or dihedral type used to measure orientation
  //         For bond, the unit vector pointing from 1 to 2 is used; for angle, the unit vector pointing from
  //         1 to 3 is used; for dihedral, the unit vector pointing from 1 to 4 is used.
  // preifx is the prefix of output files
  // noutput is the frequency to output measurement results
  // xbin,ybin,zbin are the bin numbers in x, y and z direction, respectively.
  // flag indicates how the accumulators are initialized:
  //                          0 means initialize only at first call
  //                          1 means initialize after every data output
  // qsum: if presented, the total charge on atoms in selected bond, angle, or dihedral should be this number,
  //       otherwise the bond, angle or dihedral will not be used to evaluate orientation

  if (narg < 12) error->all(FLERR,"Illegal orientation command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  orient_base = atoi(arg[4]);
  if (orient_base < 2 || orient_base > 4) error->all(FLERR, "Orient_base can only be 2 (bond), 3 (angle), or 4 (dihedral)");
  id_type = atoi(arg[5]);
  if (orient_base == 2){
    if (id_type < 1 || id_type > atom->nbondtypes) error->all(FLERR, "Wrong bond type!");
  } else if (orient_base == 3){
    if (id_type < 1 || id_type > atom->nangletypes) error->all(FLERR, "Wrong angle type!");
  } else {
    if (id_type < 1 || id_type > atom->ndihedraltypes) error->all(FLERR, "Wrong dihedral type!");
  }

  int n = strlen(arg[6])+1;
  prefix = new char[n];
  strcpy(prefix, arg[6]);

  noutput = atoi(arg[7]);

  Nbx = atoi(arg[8]);
  Nby = atoi(arg[9]);
  Nbz = atoi(arg[10]);

  restartflg = atoi(arg[11]);

  checkchg = 0;
  if (narg == 13){
     checkchg = 1;
     qref = atof(arg[12]);
  }

  first = 1;
  // enabled to output a global vector of size 6
  vector_flag = 1;
  size_vector = 6;
  global_freq = nevery;
  //extvector = 1;

}

DiagOrientation::~DiagOrientation()
{
  free(altogether);
}

int DiagOrientation::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}



void DiagOrientation::init()
{

  if (first != 1) return;
  first = 0;

#define EPSILON (1.e-10)

  subNbx=(int)ceil(((double) Nbx)/((double) comm->procgrid[0]))+12;
  subNby=(int)ceil(((double) Nby)/((double) comm->procgrid[1]))+12;
  subNbz=(int)ceil(((double) Nbz)/((double) comm->procgrid[2]))+12;
  
  sublosft[0]=domain->boxlo[0]+(double)
    ((int)((domain->sublo[0]-domain->boxlo[0])/(domain->xprd)*Nbx)-6)*
    (domain->xprd)/Nbx;
  sublosft[1]=domain->boxlo[1]+(double)
    ((int)((domain->sublo[1]-domain->boxlo[1])/(domain->yprd)*Nby)-6)*
    (domain->yprd)/Nby;
  sublosft[2]=domain->boxlo[2]+(double)
    ((int)((domain->sublo[2]-domain->boxlo[2])/(domain->zprd)*Nbz)-6)*
    (domain->zprd)/Nbz;

  subsize[0]=(domain->xprd)*(subNbx)/Nbx;
  subsize[1]=(domain->yprd)*(subNby)/Nby;
  subsize[2]=(domain->zprd)*(subNbz)/Nbz;

  //altogether = memory->create_2d_double_array(subNbx*subNby*subNbz*10+10, 10, "fix_orientation:altogether");
  altogether = (double (*)[10])calloc((size_t)(subNbx*subNby*subNbz*10+10), sizeof(double));
  if (altogether == NULL) error->one(FLERR,"Could not allocate density data\n");

  for (int i=0; i< subNbx*subNby*subNbz; i++){
    for (int j=0; j<10; j++) altogether[i][j] = 0.;
  }

  // altogether[][0] is density
  // altogether[][1-3] is momentum
  // altogether[][4-10] is second moment of momentum
  // end elements are bounds information: 3 offsets and 3 bin totals:

  if (sublosft[0] > domain->boxlo[0])
    altogether[subNbx*subNby*subNbz][4]=(double)
      ((int)((sublosft[0]+EPSILON-domain->boxlo[0])/(domain->xprd)*Nbx));
  else
    altogether[subNbx*subNby*subNbz][4]=(double)
      ((int)((sublosft[0]-EPSILON-domain->boxlo[0])/(domain->xprd)*Nbx));
  if (sublosft[1] > domain->boxlo[1])
    altogether[subNbx*subNby*subNbz][5]=(double)
      ((int)((sublosft[1]+EPSILON-domain->boxlo[1])/(domain->yprd)*Nby));
  else
    altogether[subNbx*subNby*subNbz][5]=(double)
      ((int)((sublosft[1]-EPSILON-domain->boxlo[1])/(domain->yprd)*Nby));
  if (sublosft[2] > domain->boxlo[2])
    altogether[subNbx*subNby*subNbz][6]=(double)
      ((int)((sublosft[2]+EPSILON-domain->boxlo[2])/(domain->zprd)*Nbz));
  else
    altogether[subNbx*subNby*subNbz][6]=(double)
      ((int)((sublosft[2]-EPSILON-domain->boxlo[2])/(domain->zprd)*Nbz));
  altogether[subNbx*subNby*subNbz][7]=(double)subNbx;
  altogether[subNbx*subNby*subNbz][8]=(double)subNby;
  altogether[subNbx*subNby*subNbz][9]=(double)subNbz;

  ncount = 0;
  num_ave = 0;
}


void DiagOrientation::end_of_step()
{
  double **x = atom->x;
  double  *q = atom->q;
  int nbondlist  = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;
  int nanglelist  = neighbor->nanglelist;
  int **anglelist = neighbor->anglelist;
  int ndihedrallist  = neighbor->ndihedrallist;
  int **dihedrallist = neighbor->dihedrallist;
  

  int n, i1, i2, i3, i4;
  int xb, yb, zb, place;
  double delx, dely, delz;
  double rsq, r, fwx, fwy, fwz, qsum;
  double const qzero = 0.001;
  
  cv_flag = 0;
  for (int i = 0; i < 7; ++i) vector_all[i] = vector_me[i] = 0.;

  if (orient_base == 2){ // orientation based on bond
    for (n = 0; n < nbondlist; n++) {
      if (bondlist[n][2] == id_type){
        i1 = bondlist[n][0];
        i2 = bondlist[n][1];

        if ( checkchg ){
           qsum = q[i1] + q[i2];
           if (fabs(qsum-qref) > qzero ) continue;
        }

        delx = x[i1][0] - x[i2][0];
        dely = x[i1][1] - x[i2][1];
        delz = x[i1][2] - x[i2][2];
        domain->minimum_image(delx,dely,delz);

        rsq = delx*delx + dely*dely + delz*delz;
        rsq = 1./ rsq;
        r   = sqrt(rsq);

        fwx = (x[i1][0] - delx * 0.5 - sublosft[0])/subsize[0];
        fwy = (x[i1][1] - dely * 0.5 - sublosft[1])/subsize[1];
        fwz = (x[i1][2] - delz * 0.5 - sublosft[2])/subsize[2];
        if ( fwx < 0. || fwx > 1. || fwy < 0. || fwy > 1. || fwz < 0. || fwz > 1. )
           error->warning(FLERR,"Box too small for orientation measurement!");

        xb = int(fwx*double(subNbx));
        yb = int(fwy*double(subNby));
        zb = int(fwz*double(subNbz));
        
	     place=xb*subNby*subNbz+yb*subNbz+zb;
	
        delx *= r;
        dely *= r;
        delz *= r;

        altogether[place][0] +=1.0;
        altogether[place][1] += delx*delx;
        altogether[place][2] += dely*dely;
        altogether[place][3] += delz*delz;
        altogether[place][4] += delx*dely;
        altogether[place][5] += delx*delz;
        altogether[place][6] += dely*delz;
        altogether[place][7] += delx;
        altogether[place][8] += dely;
        altogether[place][9] += delz;

        vector_me[0] +=1.0;
        vector_me[1] += delx*delx;
        vector_me[2] += dely*dely;
        vector_me[3] += delz*delz;
        vector_me[4] += delx*dely;
        vector_me[5] += delx*delz;
        vector_me[6] += dely*delz;
      }
    }
  } else if (orient_base == 3) { // orientation based on angle
    for (n = 0; n < nanglelist; n++) {
      if (anglelist[n][3] == id_type){
        i1 = anglelist[n][0];
        i2 = anglelist[n][1];
        i3 = anglelist[n][2];

        if ( checkchg ){
           qsum = q[i1] + q[i2] + q[i3];
           if ( fabs(qsum-qref) > qzero ) continue;
        }

        delx = x[i1][0] - x[i3][0];
        dely = x[i1][1] - x[i3][1];
        delz = x[i1][2] - x[i3][2];
        domain->minimum_image(delx,dely,delz);

        rsq = delx*delx + dely*dely + delz*delz;
        rsq = 1./ rsq;
        r   = sqrt(rsq);

        fwx = (x[i1][0] + delx * 0.5 - sublosft[0])/subsize[0];
        fwy = (x[i1][1] + dely * 0.5 - sublosft[1])/subsize[1];
        fwz = (x[i1][2] + delz * 0.5 - sublosft[2])/subsize[2];
        if ( fwx < 0. || fwx > 1. || fwy < 0. || fwy > 1. || fwz < 0. || fwz > 1. )
           error->warning(FLERR,"Box too small for orientation measurement!");

        xb = int(fwx*double(subNbx));
        yb = int(fwy*double(subNby));
        zb = int(fwz*double(subNbz));

        place=xb*subNby*subNbz+yb*subNbz+zb;

        delx *= r;
        dely *= r;
        delz *= r;

        altogether[place][0] +=1.;
        altogether[place][1] += delx*delx;
        altogether[place][2] += dely*dely;
        altogether[place][3] += delz*delz;
        altogether[place][4] += delx*dely;
        altogether[place][5] += delx*delz;
        altogether[place][6] += dely*delz;

        vector_me[0] +=1.0;
        vector_me[1] += delx*delx;
        vector_me[2] += dely*dely;
        vector_me[3] += delz*delz;
        vector_me[4] += delx*dely;
        vector_me[5] += delx*delz;
        vector_me[6] += dely*delz;
      }
    }
  } else if (orient_base == 4) { // orientation based on dihedral
    for (n = 0; n < ndihedrallist; n++) {
      if (dihedrallist[n][4] == id_type){
        i1 = dihedrallist[n][0];
        i2 = dihedrallist[n][1];
        i3 = dihedrallist[n][2];
        i4 = dihedrallist[n][3];
    
        if ( checkchg ){
           qsum = q[i1] + q[i2] + q[i3] + q[i4];
           if ( fabs(qsum-qref) > qzero ) continue;
        }

        delx = x[i1][0] - x[i4][0];
        dely = x[i1][1] - x[i4][1];
        delz = x[i1][2] - x[i4][2];
        domain->minimum_image(delx,dely,delz);

        rsq = delx*delx + dely*dely + delz*delz;
        rsq = 1./ rsq;
        r   = sqrt(rsq);

        fwx = (x[i1][0] + delx * 0.5 - sublosft[0])/subsize[0];
        fwy = (x[i1][1] + dely * 0.5 - sublosft[1])/subsize[1];
        fwz = (x[i1][2] + delz * 0.5 - sublosft[2])/subsize[2];
        if ( fwx < 0. || fwx > 1. || fwy < 0. || fwy > 1. || fwz < 0. || fwz > 1. )
           error->warning(FLERR,"Box too small for orientation measurement!");

        xb = int(fwx*double(subNbx));
        yb = int(fwy*double(subNby));
        zb = int(fwz*double(subNbz));

        place=xb*subNby*subNbz+yb*subNbz+zb;

        delx *= r;
        dely *= r;
        delz *= r;

        altogether[place][0] +=1.;
        altogether[place][1] += delx*delx;
        altogether[place][2] += dely*dely;
        altogether[place][3] += delz*delz;
        altogether[place][4] += delx*dely;
        altogether[place][5] += delx*delz;
        altogether[place][6] += dely*delz;
        altogether[place][7] += delx;
        altogether[place][8] += dely;
        altogether[place][9] += delz;

        vector_me[0] +=1.0;
        vector_me[1] += delx*delx;
        vector_me[2] += dely*dely;
        vector_me[3] += delz*delz;
        vector_me[4] += delx*dely;
        vector_me[5] += delx*delz;
        vector_me[6] += dely*delz;
      }
    }
  }

  num_ave++;

  if ( ++ncount == noutput ){
     ncount = 0;
    
    // get size for communication buffer
    int local_size=subNbx*subNby*subNbz*10+10, local_max, recv_size;
    int i, j, k;
    double (*buf)[10],(*dbuf)[10],(*fullset)[10], *tmp;
    MPI_Status status;
    MPI_Request request;

    MPI_Allreduce(&local_size,&local_max,1,MPI_INT,MPI_MAX,world);

    if (me==0) {
      buf = (double (*)[10]) memory->smalloc(local_max*sizeof(double),
				       "orientation:buf");

      fullset =(double (*)[10])calloc((size_t)(Nbx*Nby*Nbz*10),sizeof(double));
      if (fullset == NULL) error->one(FLERR,"Could not allocate fullset data\n");
      // in future, will set up random access temporary file to write data
      // to in place of fullset when there is not enough memory.

     for (int iproc = 0; iproc < comm->nprocs; iproc++) {
	if (iproc) {
	  MPI_Irecv(buf,local_max,MPI_DOUBLE,iproc,0,world,&request);
	  MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	  MPI_Wait(&request,&status);
	  MPI_Get_count(&status,MPI_DOUBLE,&recv_size);
	  dbuf=buf;
	} else {
	  dbuf=altogether;
	  recv_size = local_size;
	}
	recv_size /= 10;


	for (i=0; i< (int)dbuf[recv_size-1][7]; i++) {
	  xb=(int)dbuf[recv_size-1][4]+i;
	  while (xb < 0) xb+= Nbx;
	  while (xb >= Nbx) xb-= Nbx;
	  for (j=0; j< (int)dbuf[recv_size-1][8]; j++) {
	    yb=(int)dbuf[recv_size-1][5]+j;
	    while (yb < 0) yb+= Nby;
	    while (yb >= Nby) yb-= Nby;
	    for (k=0; k< (int)dbuf[recv_size-1][9]; k++) {
	      zb=(int)dbuf[recv_size-1][6]+k;
	      while (zb < 0) zb+= Nbz;
	      while (zb >= Nbz) zb-= Nbz;
	      
	      for (int idat=0; idat<10; idat++)
		fullset[(xb*Nby+yb)*Nbz+zb][idat]+=
		  dbuf[((i*(int)dbuf[recv_size-1][8]+j)*(int)dbuf[recv_size-1][9]+k)][idat];
	    }
	  }
	}
	
      }
      
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(altogether,local_size,MPI_DOUBLE,0,0,world);
    }
    
    if (me==0) {
       char outfile[128];
       FILE *outfilefp;

      sprintf(outfile,"%s."BIGINT_FORMAT,prefix,update->ntimestep);
      outfilefp= fopen(outfile,"w");
      if (outfilefp == NULL) error->one(FLERR,"Could not open density file\n");
      
      for (i=0; i<Nbx*Nby*Nbz; i++)  { 
	fprintf(outfilefp,"%16.8f ",fullset[i][0]/num_ave);
	for (int idat=1; idat<4; idat++)
	  fprintf(outfilefp,"%16.8f ",fullset[i][idat]/num_ave);
	for (int idat=4; idat<10; idat++)
	  fprintf(outfilefp,"%16.8f ",
		  fullset[i][idat]/num_ave);
	fprintf(outfilefp,"\n");
      }
      
      free(fullset);
      fclose(outfilefp);
    }
    
    // reinitialize data arrays if restartflg not zero
    if (restartflg) {
      num_ave = 0;
   	for (j=0; j< subNbx*subNby*subNbz; j++)
   	  for (int idat=0; idat<10; idat++) altogether[j][idat]=0.0;
    }
  }
  
}

double DiagOrientation::compute_vector(int n)
{
  // only sum across procs one time

  if (cv_flag == 0) {
    MPI_Allreduce(vector_me, vector_all, 7, MPI_DOUBLE, MPI_SUM, world);
    if (vector_all[0] > 0.){
      for (int i = 1; i < 7; ++i) vector_all[i] /= vector_all[0];
    }

    cv_flag = 1;
  }

  return vector_all[n+1];
}
