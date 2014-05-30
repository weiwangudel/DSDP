#if !defined(__DSDP_VMATRIXOPS_H) 
#define __DSDP_VMATRIXOPS_H
/*!
\file dsdpxmat_impl.h
\brief Structure of function pointers that each dense matrix array type 
(upper full, packed symmetric, ...) must implement.
 */

/*!
struct DSDPVMat_Ops
\brief Table of function pointers that operate on the dense matrix.
*/
struct  DSDPVMat_Ops{
  int id;
  
  int ptr_matzeroentries;
/*
../src/vecmat/dufull.c:  densematops->matzeroentries=DTRUMatZero;
../src/vecmat/dufull.c:  densematops->ptr_matzeroentries=1;
../src/vecmat/dufull.c:  densematops->matzeroentries=DTRUMatZero;
../src/vecmat/dufull.c:  densematops->ptr_matzeroentries=1;
../src/vecmat/spds.c:  dsops->matzeroentries=SpSymMatZero;
../src/vecmat/spds.c:  dsops->ptr_matzeroentries=2;
../src/vecmat/spds.c:  dsops->matzeroentries=SpSymMatZero;
../src/vecmat/spds.c:  dsops->ptr_matzeroentries=2;
../src/vecmat/diag.c:  ddiagops->matzeroentries=DiagMatZeros;
../src/vecmat/diag.c:  ddiagops->ptr_matzeroentries=3;
../src/vecmat/diag.c:  ddiagops->matzeroentries=DiagMatZeros;
../src/vecmat/diag.c:  ddiagops->ptr_matzeroentries=3;
../src/vecmat/dlpack.c:  densematops->matzeroentries=DTPUMatZero;
../src/vecmat/dlpack.c:  densematops->ptr_matzeroentries=4;
../src/vecmat/dlpack.c:  densematops->matzeroentries=DTPUMatZero;
../src/vecmat/dlpack.c:  densematops->ptr_matzeroentries=4;

*/
  
  
  int (*matgetsize)(void*,int*);
  int (*mataddouterproduct)(void*,double,double[],int);
  int (*matmult)(void*,double[],double[],int);
  int (*matscalediagonal)(void*,double);
  int (*matshiftdiagonal)(void*,double);
  int (*matfnorm2)(void*,int,double*);
  int (*matzeroentries)(void*);
  
  int (*matgeturarray)(void*,double*[],int*);
  int (*matrestoreurarray)(void*,double*[],int*);
  int (*matmineig)(void*,double[],double[],int,double*);
  int (*mattest)(void*);
  int (*matdestroy)(void*);
  int (*matview)(void*);
  const char *matname;

};

#ifdef __cplusplus
extern "C" {
#endif

extern int DSDPVMatOpsInitialize(struct  DSDPVMat_Ops*);

#ifdef __cplusplus
}
#endif

#endif


