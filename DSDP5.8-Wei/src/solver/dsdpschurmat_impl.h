#include "dsdpbasictypes.h"
#if !defined(__DSDP_SCHURMATRICES_H) 
#define __DSDP_SCHURMATRICES_H

/*!
\file dsdpschurmat_impl.h
\brief Function pointers that a Schur complement matrix (dense, sparse, parallel dense) must provide.
*/

struct  DSDPSchurMat_Ops{
  int id;
  int ptr_matzero;
  int ptr_matrownonzeros;
  /*
  ** 1: vecmat/dufull.c/DTRUMatRowNonzeros
  ** 2: vecmat/cholmat.c/DSDPGramMatRowNonzeros
  ** 3: vecmat/diag.c/DiagRowNonzeros
  ** 4: vecmat/dlpack.c/DTPUMatRowNonzeros
  */
  int ptr_mataddrow;
  int ptr_mataddelement;
  int ptr_matadddiagonal;
  int ptr_matshiftdiagonal;
  int ptr_matassemble;
  int ptr_matscaledmultiply;
  int ptr_matmultr;
  int ptr_matfactor;
  int ptr_matsolve;
  int ptr_matsetup;
  int ptr_pmatwhichdiag;
  int ptr_pmatonprocessor;
  int ptr_pmatlocalvariables;
  int ptr_pmatreduction;
  int ptr_pmatdistributed;
  int ptr_matdestroy;
  int ptr_matview;
  
  
  int (*matzero)(void*);
  int (*matrownonzeros)(void*,int,double*,int*,int);
  int (*mataddrow)(void*,int,double,double[],int);
  int (*mataddelement)(void*,int,double);
  int (*matadddiagonal)(void*,double[],int);
  int (*matshiftdiagonal)(void*,double);
  int (*matassemble)(void*);
  int (*matscaledmultiply)(void*,double[],double[],int);
  int (*matmultr)(void*,double[],double[],int);
  int (*matfactor)(void*,int*);
  int (*matsolve)(void*,double[],double[],int);
  int (*matsetup)(void*,int);
  int (*pmatwhichdiag)(void*,double[],int);
  int (*pmatonprocessor)(void*,int,int*);
  int (*pmatlocalvariables)(void*,double[],int);
  int (*pmatreduction)(void*,double[],int);
  int (*pmatdistributed)(void*,int*);
  int (*matdestroy)(void*);
  int (*matview)(void*);
  const char *matname;
};

extern int DSDPSetSchurMatOps(DSDP,struct DSDPSchurMat_Ops*,void*);
extern int DSDPSchurMatOpsInitialize(struct DSDPSchurMat_Ops*);
extern int DSDPSparsityInSchurMat(DSDP,int,int[],int);
#endif


