#if !defined(__DSDP_DATAMATRIX_H) 
#define __DSDP_DATAMATRIX_H

/*! 
\file dsdpdatamat_impl.h
\brief Structure of function pointers that each SDP data matrix type 
(sparse, dense, constant, identity, ...) must implement.
*/

/* DSDP Data Matrices have particular operations, and several implementations */
/*!
struct DSDPDataMat_Ops
\brief Table of function pointers that operate on the data matrix.
*/
struct  DSDPDataMat_Ops{
  int id;
  int ptr_matgetrank;
/*
src/vecmat/drowcol.c:  rcmatoperator->matgetrank=RCMatGetRank;
src/vecmat/drowcol.c:  rcmatoperator->ptr_matgetrank=1;
src/vecmat/vech.c:  sops->matgetrank=VechMatGetRank;
src/vecmat/vech.c:  sops->ptr_matgetrank=2;
src/vecmat/dufull.c:  sops->matgetrank=DvecumatGetRank;
src/vecmat/dufull.c:  sops->ptr_matgetrank=3;
src/vecmat/onemat.c:  cmatops->matgetrank=ConstMatGetRank;
src/vecmat/onemat.c:  cmatops->ptr_matgetrank=4;
src/vecmat/zeromat.c:  sops->matgetrank=ZGetRank;
src/vecmat/zeromat.c:  sops->ptr_matgetrank=5;
src/vecmat/identity.c:  spdiagops->matgetrank=IdentityMatGetRank;
src/vecmat/identity.c:  spdiagops->ptr_matgetrank=6;
src/vecmat/identity.c:  spdiagops->matgetrank=IdentityMatGetRank;
src/vecmat/identity.c:  spdiagops->ptr_matgetrank=6;
src/vecmat/vechu.c:  sops->matgetrank=VechMatGetRank;
src/vecmat/vechu.c:  sops->ptr_matgetrank=7;
src/vecmat/rmmat.c:  r1matops->matgetrank=R1MatGetRank;
src/vecmat/rmmat.c:  r1matops->ptr_matgetrank=8;
src/vecmat/rmmat.c:  r1matops->matgetrank=R1MatGetRank;
src/vecmat/rmmat.c:  r1matops->ptr_matgetrank=8;
src/vecmat/dlpack.c:  sops->matgetrank=DvechmatGetRank;
src/vecmat/dlpack.c:  sops->ptr_matgetrank=9;
*/
  int (*mataddallmultiple)(void*,double,double[],int,int);
  int (*matdot)(void*, double[], int, int, double *);
  int (*matgetrank)(void*,int*,int);
  int (*matgeteig)(void*,int,double*,double[],int,int[],int*);
  int (*matvecvec)(void*, double[], int,double*);
  int (*mataddrowmultiple)(void*,int,double,double[],int); /* NEEDED? */
  int (*matmultiply)(void*,double[],double[],int);
  int (*matfactor1)(void*);
  int (*matfactor2)(void*,double[],int,double[],int,double[],int,int[],int);
  int (*matfnorm2)(void*,int,double*);
  int (*matrownz)(void*,int,int[],int*,int);
  int (*matnnz)(void*,int*,int);
  int (*mattest)(void*);
  int (*matdestroy)(void*);
  int (*matview)(void*);
  const char *matname;
};

#ifdef __cplusplus
extern "C" {
#endif
extern int DSDPGetEigs(double[],int,double[],int,long int[],int, 
		       double[],int,double[],int,int[],int);
extern int DSDPGetEigs2(double[],int,double[],int,long int[],int, 
			double[],int,double[],int,int[],int);

int DSDPDataMatOpsInitialize(struct  DSDPDataMat_Ops*);

#ifdef __cplusplus
}
#endif

/*
#include "dsdpdatamat.h"
*/

#endif


