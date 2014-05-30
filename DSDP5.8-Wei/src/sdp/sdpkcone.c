#include "dsdpsdp.h"
#include "dsdpcone_impl.h"
#include "dsdpsys.h"
#include "dsdpdatamat_impl.h"   // Inlining DSDPDataMatDot Needed
//#include "dsdpdatamat.h"
#include "dsdpdualmat.h"
#include "dsdpdsmat.h"
#include "dsdpxmat.h"
#include "dsdpsys.h"
#include "dsdplanczos.h"
#include "dsdplapack.h"

#include "dsdpschurmat_impl.h"
#include "dsdpschurmat.h"
#include "dsdpbasictypes.h"

#include "dsdpxmat_impl.h"
#include "dsdpxmat.h"

#include "dsdpdualmat_impl.h"
#include "dsdpdualmat.h"

/* from vech.c and vechu.c (vechmat) */
 typedef struct {
   int    neigs;
   double *eigval;
   double *an;
   int    *cols,*nnz;
 } Eigen;
 
 typedef struct {
   int    nnzeros;
   const int    *ind;
   const double *val;
   int ishift;
   double alpha;
   
   Eigen   *Eig;
   int factored;
   int owndata;
   int    n;
 } vechmat;

#define GETI_VECH(a)    (int)(sqrt(2*(a)+0.25)-0.5)
#define GETJ_VECH(b,c)  ((b)-((c)*((c)+1))/2)

static void getij_vech(int k, int *i,int *j){
  *i=GETI_VECH(k);
  *j=GETJ_VECH(k,*i);
  return;
}

#define GETI_VECHU(a,b)    (int)((int)a/(int)b)
#define GETJ_VECHU(a,b)    (int)((int)a%(int)b)

static void getij_vechu(int k, int n, int *i,int *j){
  *i=GETI_VECHU(k,n);
  *j=GETJ_VECHU(k,n);
  return;
}

static int EigMatGetEig(Eigen* A,int row, double *eigenvalue, double eigenvector[], int n, int spind[], int *nind){
  int i,*cols=A->cols,bb,ee;
  double* an=A->an;
  *eigenvalue=A->eigval[row];
  *nind=0;
  if (cols){
    memset((void*)eigenvector,0,n*sizeof(double));
    if (row==0){ bb=0;} else {bb=A->nnz[row-1];} ee=A->nnz[row];
    for (i=bb;i<ee;i++){
      eigenvector[cols[i]]=an[i];
      spind[i-bb]=cols[i]; (*nind)++;
    }
  } else {
    memcpy((void*)eigenvector,(void*)(an+n*row),n*sizeof(double));
    for (i=0;i<n;i++)spind[i]=i;
    *nind=n;
  }
  return 0;
}

static int VechMatGetRank(void *AA,int *rank,int n){
  //printf("File %s line %d VechMatGetRank with address %d\n",__FILE__, __LINE__,&VechMatGetRank);
  vechmat*  A=(vechmat*)AA;
  switch (A->factored){
  case 1:
    *rank=A->nnzeros;
    break;
  case 2:
    *rank=2*A->nnzeros;
    break;
  case 3:
    *rank=A->Eig->neigs;
    break;
  default:
    DSDPSETERR(1,"Vech Matrix not factored yet\n");
  }
  return 0;
}

static int EigMatVecVec(Eigen* A, double v[], int n, double *vv){
  int i,rank,*cols=A->cols,neigs=A->neigs,*nnz=A->nnz,bb,ee;
  double* an=A->an,*eigval=A->eigval,dd,ddd=0;

  if (cols){
    for (rank=0;rank<neigs;rank++){
      if (rank==0){ bb=0;} else {bb=nnz[rank-1];} ee=nnz[rank];
      for (dd=0,i=bb;i<ee;i++){
	dd+=an[i]*v[cols[i]];
      }
      ddd+=dd*dd*eigval[rank];
    }
  } else {
    for (rank=0;rank<neigs;rank++){
      for (dd=0,i=0;i<n;i++){
	dd+=an[i]*v[i];
      }
      an+=n;
      ddd+=dd*dd*eigval[rank];
    }
  }
  *vv=ddd;
  return 0;
}

/* end from vech.c and vechu.c */

/* from identity.c */

typedef struct {
  int n;
  double dm;
} identitymat;

/* end from identity.c */

/* from dufull.c */
typedef enum {
  Init=0,
  Assemble=1,
  Factored=2,   /* fail to allocate required space */
  Inverted=3,    /* indefinity is detected          */
  ISymmetric=4
} MatStatus;

typedef struct{
  char UPLO;
  int LDA;
  double *val,*v2;
  double *sscale;
  double *workn;
  int  scaleit;
  int n;
  int owndata;
  MatStatus status;
} dtrumat;
/* end from dufull.c */

/* from dlpack.c */
typedef struct{
  char UPLO;
  double *val;
  double *v2;
  double *sscale;
  int  scaleit;
  int n;
  int owndata;
} dtpumat;

static void daddrow(double *v, double alpha, int i, double row[], int n){
  double *s1;
  ffinteger j,nn=n,ione=1;
  nn=i+1; s1=v+i*(i+1)/2;
  daxpy(&nn,&alpha,s1,&ione,row,&ione);
  for (j=i+1;j<n;j++){
    s1+=j;
    row[j]+=alpha*s1[i];
  }
  return;
}

/* end from dlpack.c */


/* from diag.c */
typedef struct {
  int    n;
  double *val;
  int   owndata;
} diagmat;
/* end from diag.c */
struct _P_Mat3{
  int type;
  DSDPDualMat ss;
  DSDPDSMat ds;
  SDPConeVec V;
  DSDPVMat x;
};

typedef struct _P_Mat3* Mat3;

/*! \file sdpkcone.c 
\brief Implement the DSDPCone operations using the SDPCone subroutines.
 */
static int SDPConeOperationsInitialize(struct  DSDPCone_Ops*);

static int KSDPConeSetup(void*,DSDPVec);
static int KSDPConeSetup2(void*, DSDPVec, DSDPSchurMat);
static int KSDPConeSize(void*,double*);
static int KSDPConeSparsity(void*,int,int[],int[],int);
static int KSDPConeComputeHessian(void*,double,DSDPSchurMat,DSDPVec,DSDPVec);
static int KSDPConeComputeMaxStepLength(void*, DSDPVec, DSDPDualFactorMatrix, double *);
static int KSDPConeComputeSS(void*, DSDPVec, DSDPDualFactorMatrix, DSDPTruth *);
static int KSDPConeComputeLogSDeterminant(void *, double *, double*);
static int KSDPConeComputeXX(void*, double, DSDPVec,DSDPVec,DSDPVec,double*);
static int KSDPConeDestroy(void*);

#undef __FUNCT__  
#define __FUNCT__ "KSDPConeComputeHessian"
static int KSDPConeComputeHessian( void *K, double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2){
  //printf("sdp/sdpkcone.c/KSDPConeComputeHessian with address %d\n",&KSDPConeComputeHessian);
  int info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  //info=SDPConeComputeHessian(sdpcone,mu,M,vrhs1,vrhs2);DSDPCHKERR(info);
  //Wei: inline (2) Function Body of SDPConeComputeHessian
  // Parameter names are the same, so just wrap all things inside the new {}
  {
  //SDPConeComputeHessian( SDPCone sdpcone, double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2){
    int i,k,kt,kk,m,n,rank,info;
    int ncols,ii,idA;
    double rtemp,ack,ggamma,bmu,scl;
    double rhs1i,rhs2i;
    DSDPTruth method1;
    SDPConeVec W,W2;
    DSDPVec MRowI=sdpcone->Work, Select=sdpcone->Work2;
    DSDPDataMat AA;
    DSDPDualMat S;
    DSDPVMat T;
    DSDPDataTranspose ATranspose=sdpcone->ATR;
    SDPblk *blk=sdpcone->blk;
    DSDPIndex IS;
  
    /* Evaluate M */
   //Wei  DSDPFunctionBegin;
   //Wei  SDPConeValid(sdpcone);
   //Wei  info=DSDPVecGetSize(vrhs1,&m);DSDPCHKERR(info);
    m = vrhs1.dim;
  
    for (i=0; i<m; i++){  /* One row at a time */
  
      /* Which Coluns */
      rhs1i=0;rhs2i=0;
      //info=DSDPVecZero(MRowI);DSDPCHKERR(info);
      {
      //int DSDPVecZero(DSDPVec V){
	DSDPVec V=MRowI;
        int n=V.dim;
        double *v=V.val;
        memset((void*)v,0,n*sizeof(double));
      //}

      } // end of DSDPVecZero
      //info=DSDPSchurMatRowColumnScaling(M,i,Select,&ncols); DSDPCHKERR(info); 
      {
      //int DSDPSchurMatRowColumnScaling(DSDPSchurMat M,int row, DSDPVec V, int *nzcols){
	int row = i;
        DSDPVec V=Select;
        int *nzcols = &ncols;

        int info,m;
        double *cols,r=M.schur->r;
        DSDPTruth flag;
        //info=DSDPVecSet(0.0,V);DSDPCHKERR(info);
	{
        //int DSDPVecSet(double alpha, DSDPVec V){
       	  double alpha = 0.0; 
          int i,ii,n=V.dim;
          double *val=V.val;
        
          if (alpha==0.0){
            memset((void*)val,0,n*sizeof(double));
          }
          /* for (ii=0; ii<n/4; ++ii){
            i=ii*4;
            val[i] = val[i+1] = val[i+2] = val[i+3] = alpha; 
          }
          for (i=4*(n/4); i<n; ++i){
            val[i]= alpha;
          } */
        //}
	}  // end of DSDPVecSet
        info=DSDPVecGetSize(V,&m);DSDPCHKERR(info);
        if (row==0){
	  //info=DSDPVecZero(V);DSDPCHKERR(info);
	  {
            int n=V.dim;
            double *v=V.val;
            memset((void*)v,0,n*sizeof(double));
	  } // end of DSDPVecZero
	  *nzcols=0;
	}
        else if (row==m-1){
          //info=DSDPVecZero(V);DSDPCHKERR(info);
	  {
  	    int n=V.dim;
  	    double *v=V.val;
  	    memset((void*)v,0,n*sizeof(double));
	  } // end of DSDPVecZero
	  *nzcols=0;
          if (r){info=DSDPVecSetR(V,1.0);DSDPCHKERR(info);*nzcols=1;}
        } else if (M.dsdpops->matrownonzeros){
          //printf("M.dsdpops->ptr_matrownonzeros %d\n",M.dsdpops->ptr_matrownonzeros);
          if (M.dsdpops->ptr_matrownonzeros != 1) {
	    printf("Not a static pointer functin, Do not just map to DTRUMatRowNonzeros!!!!!!!\n");
	    exit(-1);
	  }
          info=DSDPVecGetSize(V,&m);DSDPCHKERR(info);
          info=DSDPVecGetArray(V,&cols);DSDPCHKERR(info);
          //printf("File %s line %d M.dsdpops->matrownonzeros point to %d located in ",__FILE__, __LINE__,M.dsdpops->matrownonzeros);
          //Wei commented!!!: info=(M.dsdpops->matrownonzeros)(M.data,row-1,cols+1,nzcols,m-2); //DSDPChkMatError(M,info);
          

          //info=(M.dsdpops->matrownonzeros)(M.data,row-1,cols+1,nzcols,m-2); //DSDPChkMatError(M,info);
          {
          //static int DTRUMatRowNonzeros(void*M, int row, double cols[], int *ncols,int nrows){
            int i;
            *nzcols = (row-1)+1;
            for (i=0;i<=(row-1);i++){
              (cols+1)[i]=1.0;
            }
            memset((void*)((cols+1)+(row-1)+1),0,(m-2-(row-1)-1)*sizeof(int));
          //}
	  } // end inlining function pointer matrownonzeros          


          info=DSDPVecRestoreArray(V,&cols);DSDPCHKERR(info);
          //info=DSDPZeroFixedVariables(M,V);DSDPCHKERR(info);
	  {
          //int DSDPZeroFixedVariables( DSDPSchurMat M, DSDPVec dy){
            int i,info; 
            FixedVariables *fv=&M.schur->fv;
            for (i=0;i<fv->nvars;i++){
              info=DSDPVecSetElement(V,fv->var[i],0.0);DSDPCHKERR(info);
            }
          //}


	  }  // end of DSDPZeroFixedVariables 
          info=DSDPVecSetC(V,0.0);DSDPCHKERR(info);
          if (r){info=DSDPVecSetR(V,1.0);DSDPCHKERR(info);}
          //info=DSDPIsFixed(M,row,&flag);DSDPCHKERR(info); 
	  {
	  //int DSDPIsFixed( DSDPSchurMat M, int vari, DSDPTruth *flag){
	    int vari = row;
	    int i;
	    FixedVariables *fv=&M.schur->fv;
	    flag=DSDP_FALSE;
	    for (i=0;i<fv->nvars;i++){
	      if (fv->var[i]==vari){
	        flag=DSDP_TRUE;
	        break;
	      }
	    }
	  //} original end 

	  } // end of DSDPIsFixed 
		
          if (flag==DSDP_TRUE&&*nzcols>0){info=DSDPVecZero(V);*nzcols=0;DSDPFunctionReturn(0);}
        } else {
          //DSDPNoOperationError(M);
	  exit(-1);
        }
      //}
	
      } // end of DSDPSchurMatRowColumnScaling

      if (ncols==0){continue;}
   
      for (kt=0; kt<ATranspose.nnzblocks[i]; kt++){ /* Loop over all blocks */
        kk=ATranspose.nzblocks[i][kt];
        idA=ATranspose.idA[i][kt];
        //info=DSDPBlockGetMatrix(&blk[kk].ADATA,idA,&ii,&scl,&AA);DSDPCHKBLOCKERR(kk,info);
	{
        //int DSDPBlockGetMatrix(DSDPBlockData *ADATA,int id, int *vari, double *scl, DSDPDataMat *A){
	  DSDPBlockData *ADATA = &blk[kk].ADATA;
          int id = idA;
	  int *vari = &ii;
	  DSDPDataMat *A=&AA;

          if (id>=0 && id < ADATA->nnzmats){
            if (vari) *vari=ADATA->nzmat[id];
            if (scl) scl=ADATA->scl;
            if (A) *A=ADATA->A[id];
          } else {
            DSDPSETERR2(2,"Invalid Matrix request.  0 <= %d < %d\n",id,ADATA->nnzmats);
          }
          //}
	} // end of DSDPBlockGetMatrix 
        //Wei      if (ii!=i){DSDPSETERR2(8,"Data Transpose Error: var %d does not equal %d.\n",i,ii);}
        //info = DSDPDataMatGetRank(AA,&rank,blk[kk].n);DSDPCHKBLOCKERR(kk,info);
	{
        //int DSDPDataMatGetRank(DSDPDataMat A, int *rank, int n){
	  DSDPDataMat A=AA;
	  int n = blk[kk].n; 
          int info;
          if (A.dsdpops->matgetrank){
            //printf("A.dsdpops->matgetrank %d\n",A.dsdpops->ptr_matgetrank);
            //printf("File %s line %d A.dsdpops->matgetrank point to %d located in ",__FILE__, __LINE__,A.dsdpops->matgetrank);
            //src/vecmat/vech.c:  sops->matgetrank=VechMatGetRank;
            //src/vecmat/vech.c:  sops->ptr_matgetrank=2;
            //src/vecmat/identity.c:  spdiagops->matgetrank=IdentityMatGetRank;
            //src/vecmat/identity.c:  spdiagops->ptr_matgetrank=6;
            //src/vecmat/identity.c:  spdiagops->matgetrank=IdentityMatGetRank;
            //src/vecmat/identity.c:  spdiagops->ptr_matgetrank=6; 
            //src/vecmat/vechu.c:  sops->matgetrank=VechMatGetRank;
            //src/vecmat/vechu.c:  sops->ptr_matgetrank=7;
            
            if (A.dsdpops->ptr_matgetrank == 2 || A.dsdpops->ptr_matgetrank == 7 ) { /* vech.c and vechu.c the same) */
            //info=(A.dsdpops->matgetrank)(A.matdata,&rank,n); //DSDPChkDataError(A,info);
	      { 
              //static int VechMatGetRank(void *AA,int *rank,int n)
                //printf("File %s line %d VechMatGetRank with address %d\n",__FILE__, __LINE__,&VechMatGetRank);
                vechmat*  B=(vechmat*)A.matdata;
                switch (B->factored){
                case 1:
                  rank=B->nnzeros;
                  break;
                case 2:
                  rank=2*B->nnzeros;                                               
                  break;                
                case 3:
                  rank=B->Eig->neigs;
                  break;
                default:
                  DSDPSETERR(1,"Vech Matrix not factored yet\n");
                } // end of switch
              }   // end of inlining VechMatGetRank
	   
            } 
            else if (A.dsdpops->ptr_matgetrank == 6) {
              //info=(A.dsdpops->matgetrank)(A.matdata,&rank,n); //DSDPChkDataError(A,info);
            {
              //static int IdentityMatGetRank(void *AA, int*rank, int n){
                  //printf("File %s line %d IdentityMatGetRank with address %d\n",__FILE__, __LINE__,&IdentityMatGetRank);
                  identitymat* B=(identitymat*)A.matdata;
                  rank=B->n;
                  
            }  
	    }
            else {
	      printf(" Error! Should be impossible to get here!\n");
	      exit(-1);
	    } 
            
            //info=(A.dsdpops->matgetrank)(A.matdata,&rank,n); //DSDPChkDataError(A,info);
            
          } else {
            //DSDPNoOperationError(A);
	    exit(-1);
          }
        //} // original  end
	} // end of DSDPDataMatGetRank
        if (rank==0) continue;
  
        T=blk[kk].T; S=blk[kk].S; W=blk[kk].W; W2=blk[kk].W2;
        n=blk[kk].n; ggamma=blk[kk].gammamu; bmu=blk[kk].bmu; IS=blk[kk].IS;
  
        method1=DSDP_TRUE;  /* Simple heuristic */
        if (rank==1) method1=DSDP_FALSE;
        if (rank==2 && ncols<=n) method1=DSDP_FALSE;
        if (rank*rank*ncols<=n+1)method1=DSDP_FALSE;
        if (ncols*blk[kk].nnz < n*n/10) method1=DSDP_FALSE;
        if (ncols==1 && i==m-1)method1=DSDP_FALSE;
        if (n<5) method1=DSDP_TRUE;
       //wei      if (0==1) method1=DSDP_FALSE;
        if (method1==DSDP_TRUE){
	  //info=DSDPVMatZeroEntries(T);DSDPCHKBLOCKERR(kk,info);
          {
          //int DSDPVMatZeroEntries(DSDPVMat X){
	    DSDPVMat X=T; 
            int info;
            if (X.dsdpops->matzeroentries){
            //info=(X.dsdpops->matzeroentries)(X.matdata); //DSDPChkMatError(X,info);
            //printf("X.dsdpops->matzeroentries %d\n",X.dsdpops->ptr_matzeroentries);
            // only 1 and 4
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
            if(X.dsdpops->ptr_matzeroentries==1){
            //info=(X.dsdpops->matzeroentries)(X.matdata); //DSDPChkMatError(X,info);
            //static int DTRUMatZero(void* AA){
            {
              //printf("File %s line %d DTRUMatZero with address %d\n",__FILE__, __LINE__,&DTRUMatZero);
              dtrumat* B=(dtrumat*) X.matdata;
              int mn=B->n*(B->LDA);
              double *vv=B->val;
              memset((void*)vv,0,mn*sizeof(double));
              B->status=Assemble;
              //return 0;
            }
            }else if(X.dsdpops->ptr_matzeroentries==4){
              //printf("File %s line %d X.dsdpops->matzeroentries point to %d located in ",__FILE__, __LINE__,X.dsdpops->matzeroentries);
              //info=(X.dsdpops->matzeroentries)(X.matdata); //DSDPChkMatError(X,info);
            //static int DTPUMatZero(void* AA){
            {
              //printf("File %s line %d DTPUMatZero with address %d\n",__FILE__, __LINE__,&DTPUMatZero);
              dtpumat* B=(dtpumat*) X.matdata;
              int mn=B->n*(B->n+1)/2;
              double *vv=B->val;
              memset((void*)vv,0,mn*sizeof(double));
              //return 0;
            }
            }else{
                printf("X.dsdpops->ptr_matzeroentries Error! Should be impossible to get here!\n");
            }
              
              
            } else {
	      exit(-1);
              //DSDPNoOperationError(X);
            }
          //}
	  } // end of inlining DSDPVMatZeroEntries
	}
        for (k=0; k<rank; k++){
  	
  	//info=DSDPDataMatGetEig(AA,k,W,IS,&ack); DSDPCHKBLOCKERR(kk,info);
	{
        //int DSDPDataMatGetEig(DSDPDataMat A, int rr, SDPConeVec V, DSDPIndex S, double *eigenvalue){
	  DSDPDataMat A=AA;
	  int rr = k;
	  SDPConeVec V=W;
	  DSDPIndex S=IS;
	  double *eigenvalue =&ack;

          int info,n;
          double *vv;
          if (A.dsdpops->matgeteig){
            //printf("A.dsdpops->matgeteig %d\n",A.dsdpops->ptr_matgeteig);
            // only 2,6 and 7
            
            /*
            ../src/vecmat/drowcol.c:  rcmatoperator->matgeteig=RCMatGetEig;
            ../src/vecmat/drowcol.c:  rcmatoperator->ptr_matgeteig=1;
            ../src/vecmat/vech.c:  sops->matgeteig=VechMatGetEig;
            ../src/vecmat/vech.c:  sops->ptr_matgeteig=2;
            ../src/vecmat/dufull.c:  sops->matgeteig=DvecumatGetEig;
            ../src/vecmat/dufull.c:  sops->ptr_matgeteig=3;
            ../src/vecmat/onemat.c:  cmatops->matgeteig=ConstMatGetEig;
            ../src/vecmat/onemat.c:  cmatops->ptr_matgeteig=4;
            ../src/vecmat/zeromat.c:  sops->matgeteig=ZGetEig;
            ../src/vecmat/zeromat.c:  sops->ptr_matgeteig=5;
            ../src/vecmat/identity.c:  spdiagops->matgeteig=IdentityMatGetEig;
            ../src/vecmat/identity.c:  spdiagops->ptr_matgeteig=6;
            ../src/vecmat/identity.c:  spdiagops->matgeteig=IdentityMatGetEig;
            ../src/vecmat/identity.c:  spdiagops->ptr_matgeteig=6;
            ../src/vecmat/vechu.c:  sops->matgeteig=VechMatGetEig;
            ../src/vecmat/vechu.c:  sops->ptr_matgeteig=7;
            ../src/vecmat/rmmat.c:  r1matops->matgeteig=R1MatGetEig;
            ../src/vecmat/rmmat.c:  r1matops->ptr_matgeteig=8;
            ../src/vecmat/rmmat.c:  r1matops->matgeteig=R1MatGetEig;
            ../src/vecmat/rmmat.c:  r1matops->ptr_matgeteig=8;
            ../src/vecmat/dlpack.c:  sops->matgeteig=DvechmatGetEig;
            ../src/vecmat/dlpack.c:  sops->ptr_matgeteig=9;
            */ 
            
            
            
            //printf("File %s line %d A.dsdpops->matgeteig point to %d located in ",__FILE__, __LINE__,A.dsdpops->matgeteig);
            info=SDPConeVecGetArray(V,&vv); DSDPCHKERR(info);
            info=SDPConeVecGetSize(V,&n); DSDPCHKERR(info);
            //info=(A.dsdpops->matgeteig)(A.matdata,rr, eigenvalue, vv,n,S.indx+1,S.indx); //DSDPChkDataError(A,info);
            if(A.dsdpops->ptr_matgeteig == 2){
            {
              //info=(A.dsdpops->matgeteig)(A.matdata,rr, eigenvalue, vv,n,S.indx+1,S.indx); //DSDPChkDataError(A,info);
            //static int VechMatGetEig(void* AA, int rank, double *eigenvalue, double vv[], int n, int indx[], int *nind){
              //printf("File %s line %d VechMatGetEig with address %d\n",__FILE__, __LINE__,&VechMatGetEig);
              vechmat*  B=(vechmat*)A.matdata;
              const double *val=B->val,tt=sqrt(0.5);
              int info,i,j,k,ishift=B->ishift;
              const int *ind=B->ind;

              S.indx[0]=0;
              switch (B->factored){
              case 1:
                memset(vv,0,n*sizeof(double));
                getij_vech(ind[rr]-ishift,&i,&j);
                vv[i]=1.0;
                *eigenvalue=val[rr]*B->alpha;
                S.indx[0]=1;
                S.indx[1]=i;
                break;
              case 2:
                memset(vv,0,n*sizeof(double));
                k=rr/2;
                getij_vech(ind[k]-ishift,&i,&j);
                if (i==j){
                  if (k*2==rr){ 
	            vv[i]=1.0; *eigenvalue=val[k]*B->alpha;
	            S.indx[0]=1;
	            S.indx[1]=i;
                  } else {
	            *eigenvalue=0;
                  }
                } else {
                  if (k*2==rr){
	            vv[i]=tt;  vv[j]=tt; *eigenvalue=val[k]*B->alpha;
	            S.indx[0]=2;
	            S.indx[1]=i; S.indx[2]=j;
                  } else {
	            vv[i]=-tt; vv[j]=tt; *eigenvalue=-val[k]*B->alpha;
	            S.indx[0]=2;
	            S.indx[1]=i; S.indx[2]=j;
                  }
                }
                break;
              case 3:
                info=EigMatGetEig(B->Eig,rr,eigenvalue,vv,n,S.indx+1,S.indx);DSDPCHKERR(info);
                //EigMatGetEig(Eigen* A,int row, double *eigenvalue, double eigenvector[], int n, int spind[], int *nind)
                  
                *eigenvalue=*eigenvalue*B->alpha;
                break;
              default:
                DSDPSETERR(1,"Vech Matrix not factored yet\n");
              }

              //return 0;  
            }
            }else if(A.dsdpops->ptr_matgeteig == 6){
            
            {
            //info=(A.dsdpops->matgeteig)(A.matdata,rr, eigenvalue, vv,n,S.indx+1,S.indx); //DSDPChkDataError(A,info);
            //static int IdentityMatGetEig(void*AA, int neig, double *eig, double v[], int n, int* indx, int *nind){
              //printf("File %s line %d IdentityMatGetEig with address %d\n",__FILE__, __LINE__,&IdentityMatGetEig);
              identitymat* B = (identitymat*)A.matdata;

              if (rr<0 || rr>= B->n){ *eigenvalue=0; return 0;} 
              memset((void*)vv,0,(B->n)*sizeof(double)); 
              vv[rr]=1.0;
              S.indx[1]=rr;
              S.indx[0]=1;
              *eigenvalue=B->dm;  
              //return 0;
            }
            
            
            }else if(A.dsdpops->ptr_matgeteig == 7){
                //info=(A.dsdpops->matgeteig)(A.matdata,rr, eigenvalue, vv,n,S.indx+1,S.indx); //DSDPChkDataError(A,info);
                //static int VechMatGetEig(void* AA, int rank, double *eigenvalue, double vv[], int n, int indx[], int *nind){
                {  //printf("File %s line %d VechMatGetEig with address %d\n",__FILE__, __LINE__,&VechMatGetEig);
                  vechmat*  B=(vechmat*)A.matdata;
                  const double *val=B->val,tt=sqrt(0.5);
                  int info,i,j,k,t;
                  const int *ind=B->ind,ishift=B->ishift;

                  S.indx[0]=0;
                  switch (B->factored){
                  case 1:
                    memset(vv,0,n*sizeof(double));
                    t=ind[rr]-ishift;
                    i=GETI_VECHU(t,n);
                    j=GETJ_VECHU(t,n);
                    vv[i]=1.0;
                    *eigenvalue=val[rr]*B->alpha;
                    S.indx[0]=1;
                    S.indx[1]=i;
                    break;
                  case 2:
                    memset(vv,0,n*sizeof(double));
                    k=rr/2;
                    getij_vechu(ind[k]-ishift,n,&i,&j);
                    if (i==j){
                      if (k*2==rr){ 
	                vv[i]=1.0; *eigenvalue=val[k]*B->alpha;
	                S.indx[0]=1;
	                S.indx[1]=i;
                      } else {
	                *eigenvalue=0;
                      }
                    } else {
                      if (k*2==rr){
	                vv[i]=tt;  vv[j]=tt; *eigenvalue=val[k]*B->alpha;
	                S.indx[0]=2;
	                S.indx[1]=i; S.indx[2]=j;
                      } else {
	                vv[i]=-tt; vv[j]=tt; *eigenvalue=-val[k]*B->alpha;
	                S.indx[0]=2;
	                S.indx[1]=i; S.indx[2]=j;
                      }
                    }
                    break;
                  case 3:
                    info=EigMatGetEig(B->Eig,rr,eigenvalue,vv,n,S.indx+1,S.indx);DSDPCHKERR(info);
                    *eigenvalue=*eigenvalue*B->alpha;
                    break;
                  default:
                    DSDPSETERR(1,"Vech Matrix not factored yet\n");
                  }

                  //return 0;  
                }
            }
            
            
            else{
                printf("A.dsdpops->ptr_matgeteig Error! Should be impossible to get here!\n");
                //info=(A.dsdpops->matgeteig)(A.matdata,rr, eigenvalue, vv,n,S.indx+1,S.indx); //DSDPChkDataError(A,info);
            }
            
            
            info=SDPConeVecRestoreArray(V,&vv); DSDPCHKERR(info);
          } else {
            //DSDPNoOperationError(A);
	    exit (-1);
          }
        //}


	} // end of DSDPDataMatGetEig
  	if (ack==0.0) continue;
  	ack*=scl;
  	//info=DSDPDualMatInverseMultiply(S,IS,W,W2);DSDPCHKBLOCKERR(kk,info);
	{
        //int DSDPDualMatInverseMultiply(DSDPDualMat S, DSDPIndex IS, SDPConeVec B, SDPConeVec X){
	  SDPConeVec B=W;
          SDPConeVec X=W2;
          int info,n;
          double *bb,*xx;
          //DSDPEventLogBegin(sdpdualsolve);
          if (S.dsdpops->matinversemultiply){
            
            //printf("File %s line %d S.dsdpops->matinversemultiply point to %d located in ",__FILE__, __LINE__,S.dsdpops->matinversemultiply);
            info=SDPConeVecGetSize(X,&n); DSDPCHKERR(info);
            info=SDPConeVecGetArray(B,&bb); DSDPCHKERR(info);
            info=SDPConeVecGetArray(X,&xx); DSDPCHKERR(info);
            //info=(S.dsdpops->matinversemultiply)(S.matdata,IS.indx+1,IS.indx[0],bb,xx,n); //DSDPChkDMatError(S,info);
            //printf("S.dsdpops->matinversemultiply %d\n",S.dsdpops->ptr_matinversemultiply);
            //only 3 and 4
            /*
            ../src/vecmat/dufull.c:  sops->matinversemultiply=DTRUMatInverseMultiply;
            ../src/vecmat/dufull.c:  sops->ptr_matinversemultiply=1;
            ../src/vecmat/dufull.c:  sops->matinversemultiply=DTRUMatInverseMultiply;
            ../src/vecmat/dufull.c:  sops->ptr_matinversemultiply=1;
            ../src/vecmat/cholmat2.c:  sops->matinversemultiply=SMatSolve;
            ../src/vecmat/cholmat2.c:  sops->ptr_matinversemultiply=2;
            ../src/vecmat/diag.c:  sops->matinversemultiply=DiagMatSolve2;
            ../src/vecmat/diag.c:  sops->ptr_matinversemultiply=3;
            ../src/vecmat/diag.c:  sops->matinversemultiply=DiagMatSolve2;
            ../src/vecmat/diag.c:  sops->ptr_matinversemultiply=3;
            ../src/vecmat/dlpack.c:  sops->matinversemultiply=DTPUMatInverseMult;
            ../src/vecmat/dlpack.c:  sops->ptr_matinversemultiply=4;
            */
            if(S.dsdpops->ptr_matinversemultiply==3){
                //info=(S.dsdpops->matinversemultiply)(S.matdata,IS.indx+1,IS.indx[0],bb,xx,n); //DSDPChkDMatError(S,info);
                //static int DiagMatSolve2(void* A, int indx[], int nindx, double b[], double x[],int n){
                {  
                  //printf("File %s line %d DiagMatSolve with address %d\n",__FILE__, __LINE__,&DiagMatSolve);
                  diagmat* AA = (diagmat*)S.matdata;
                  double *v=AA->val;
                  int i,j;
                  memset((void*)xx,0,n*sizeof(double));
                  for (j=0;j<IS.indx[0];j++){
                    i=IS.indx[1+j];
                    xx[i]=bb[i]/v[i];
                  }
                  //return 0;
                }
                
                
            
            }else if(S.dsdpops->ptr_matinversemultiply==4){
                info=(S.dsdpops->matinversemultiply)(S.matdata,IS.indx+1,IS.indx[0],bb,xx,n); //DSDPChkDMatError(S,info);
                //static int DTPUMatInverseMult(void* AA, int indx[], int nind, double x[], double y[], int n){
                {  
                  //printf("File %s line %d DTPUMatInverseMult with address %d\n",__FILE__, __LINE__,&DTPUMatInverseMult);
                  dtpumat* A=(dtpumat*) S.matdata;
                  ffinteger ione=1,N=n;
                  double BETA=0.0,ALPHA=1.0;
                  double *AP=A->v2,*Y=xx,*X=bb;
                  int i,ii;
                  char UPLO=A->UPLO;
                  
                  if (A->n != n) info= 1;
                  if (bb==0 && n>0) info= 3;
                  
                  if (IS.indx[0]<n/4 ){
                    memset((void*)xx,0,n*sizeof(double));    
                    for (ii=0;ii<IS.indx[0];ii++){
                      i=IS.indx[1+ii];  ALPHA=bb[i];
                      daddrow(AP,ALPHA,i,xx,n);
                    }
                  } else {
                    ALPHA=1.0;
                    dspmv(&UPLO,&N,&ALPHA,AP,X,&ione,&BETA,Y,&ione);
                  }
                  info= 0;
                }
            
            }else{
                printf("S.dsdpops->ptr_matinversemultiply Error! Should be impossible to get here!\n");
            }
            
            
            
            info=SDPConeVecRestoreArray(X,&xx); DSDPCHKERR(info);
            info=SDPConeVecRestoreArray(B,&bb); DSDPCHKERR(info);
          } else {
            //DSDPNoOperationError(S);
	    exit (-1);
          }
          //DSDPEventLogEnd(sdpdualsolve);
        //}  end of original

	} // end of DSDPDualMatInverseMultiply
  
  	/* RHS terms */
  	info = SDPConeVecDot(W,W2,&rtemp); DSDPCHKBLOCKERR(kk,info);
	//{
        ////int SDPConeVecDot(SDPConeVec V1, SDPConeVec V2, double *ans){
	//  double *ans = &rtemp;
	//  SDPConeVec V1 = W;
	//  SDPConeVec V2 = W2;
        //  ffinteger ione=1, nn=V1.dim;
        //  double *v1=V1.val,*v2=V2.val;
        //  *ans=ddot(&nn,v1,&ione,v2,&ione);
        ////}
	//} // end of SDPConeVecDot
  	if (rtemp==0.0) continue;
  	rhs1i+=rtemp*ack*bmu; rhs2i+=rtemp*ack*ggamma*mu;
  	ack*=(ggamma+bmu);
  
  	if (method1==DSDP_TRUE){
  	  //info=DSDPVMatAddOuterProduct(T,ack*mu,W2);DSDPCHKBLOCKERR(kk,info);
	  {
          //int DSDPVMatAddOuterProduct(DSDPVMat X, double alpha, SDPConeVec V){
	    DSDPVMat X=T;
            double alpha=ack*mu;
	    SDPConeVec V=W2;

            int info,n;
            double *v;
            //DSDPEventLogBegin(sdpxmatevent);
            info=SDPConeVecGetSize(V,&n); //DSDPCHKERR(info);
            if (X.dsdpops->mataddouterproduct){
              //printf("X.dsdpops->mataddouterproduct %d\n",X.dsdpops->ptr_mataddouterproduct);
            
              //printf("File %s line %d X.dsdpops->mataddouterproduct point to %d located in ",__FILE__, __LINE__,X.dsdpops->mataddouterproduct);
              info=SDPConeVecGetArray(V,&v); DSDPCHKERR(info);
              
              //info=(X.dsdpops->mataddouterproduct)(X.matdata,alpha,v,n); //DSDPChkMatError(X,info);
              /*../src/vecmat/dufull.c:  densematops->mataddouterproduct=DTRUMatOuterProduct;
                ../src/vecmat/dufull.c:  densematops->ptr_mataddouterproduct=1;
                ../src/vecmat/dlpack.c:  densematops->mataddouterproduct=DTPUMatOuterProduct;
                ../src/vecmat/dlpack.c:  densematops->ptr_mataddouterproduct=2;*/
                
              if(X.dsdpops->ptr_mataddouterproduct == 1){
                //info=(X.dsdpops->mataddouterproduct)(X.matdata,alpha,v,n); //DSDPChkMatError(X,info);
                //static int DTRUMatOuterProduct(void* AA, double alpha, double x[], int n){
                {  
                  //printf("File %s line %d DTRUMatOuterProduct with address %d\n",__FILE__, __LINE__,&DTRUMatOuterProduct);
                  dtrumat* B=(dtrumat*) X.matdata;
                  ffinteger ione=1,N=n,LDA=B->LDA;
                  double *vv=B->val;
                  char UPLO=B->UPLO;
                  dsyr(&UPLO,&N,&alpha,v,&ione,vv,&LDA);
                  //return 0;
                }
              
              }else if(X.dsdpops->ptr_mataddouterproduct == 2){
                //info=(X.dsdpops->mataddouterproduct)(X.matdata,alpha,v,n); //DSDPChkMatError(X,info);
                //static int DTPUMatOuterProduct(void* AA, double alpha, double x[], int n){
                {  
                  //printf("File %s line %d DTPUMatOuterProduct with address %d\n",__FILE__, __LINE__,&DTPUMatOuterProduct);
                  dtpumat* B=(dtpumat*) X.matdata;
                  ffinteger ione=1,N=n;
                  double *vv=B->val;
                  char UPLO=B->UPLO;
                  dspr(&UPLO,&N,&alpha,v,&ione,vv);
                  //return 0;
                }
              
              
              }else{
                printf("X.dsdpops->ptr_mataddouterproduct Error!");
                
              }
              
              
              info=SDPConeVecRestoreArray(V,&v); //DSDPCHKERR(info);
            } else {
              //DSDPNoOperationError(X);
	      exit (-1);
            }
            //DSDPEventLogEnd(sdpxmatevent);
          //}
	  } // end of DSDPVMatAddOuterProduct
  	} else {
  	  //info=DSDPBlockvAv(&blk[kk].ADATA,ack*mu,Select,W2,MRowI);//DSDPCHKBLOCKERR(kk,info);
	  {
          //int DSDPBlockvAv(DSDPBlockData *ADATA, double aa, DSDPVec Alpha, SDPConeVec V, DSDPVec VAV){
          
            DSDPBlockData *ADATA = &blk[kk].ADATA;
            int    ii,vari,info;
            double sum=0,aalpha=0,scl=ADATA->scl;
            double aa = ack*mu;
            DSDPVec Alpha = Select;
            SDPConeVec V = W2;
            DSDPVec VAV = MRowI;            
 
            //DSDPEventLogBegin(sdpvecvecevent);
            if (aa==0){/* do nothing */ }
	    else {
            for (ii=0;ii<ADATA->nnzmats; ii++){  /* Matrix Entries */
              vari=ADATA->nzmat[ii];
              info=DSDPVecGetElement(Alpha,vari,&aalpha);DSDPCHKVARERR(vari,info);
              if (aalpha==0.0) continue;
              //info=DSDPDataMatVecVec(ADATA->A[ii],V,&sum);DSDPCHKVARERR(vari,info);
	      {
              //int DSDPDataMatVecVec(DSDPDataMat A, SDPConeVec W, double *v){
		DSDPDataMat A = ADATA->A[ii];
		SDPConeVec W = V;
		double*v = &sum;
                int info,n;
                double *x;
              
                if (A.dsdpops->matvecvec){
                //printf("A.dsdpops->matvecvec %d\n",A.dsdpops->ptr_matvecvec);
                // only 2 and 6
                /*
                ../src/vecmat/drowcol.c:  rcmatoperator->matvecvec=RCMatVecVec;
                ../src/vecmat/drowcol.c:  rcmatoperator->ptr_matvecvec=1;
                ../src/vecmat/vech.c:  sops->matvecvec=VechMatVecVec;
                ../src/vecmat/vech.c:  sops->ptr_matvecvec=2;
                ../src/vecmat/dufull.c:  sops->matvecvec=DvecumatVecVec;
                ../src/vecmat/dufull.c:  sops->ptr_matvecvec=3;
                ../src/vecmat/onemat.c:  cmatops->matvecvec=ConstMatVecVec;
                ../src/vecmat/onemat.c:  cmatops->ptr_matvecvec=4;
                ../src/vecmat/zeromat.c:  sops->matvecvec=ZVecVec;
                ../src/vecmat/zeromat.c:  sops->ptr_matvecvec=5;
                ../src/vecmat/identity.c:  spdiagops->matvecvec=IdentityMatVecVec;
                ../src/vecmat/identity.c:  spdiagops->ptr_matvecvec=6;
                ../src/vecmat/identity.c:  spdiagops->matvecvec=IdentityMatVecVec;
                ../src/vecmat/identity.c:  spdiagops->ptr_matvecvec=6;
                ../src/vecmat/vechu.c:  sops->matvecvec=VechMatVecVec;
                ../src/vecmat/vechu.c:  sops->ptr_matvecvec=7;
                ../src/vecmat/rmmat.c:  r1matops->matvecvec=R1MatVecVec;
                ../src/vecmat/rmmat.c:  r1matops->ptr_matvecvec=8;
                ../src/vecmat/rmmat.c:  r1matops->matvecvec=R1MatVecVec;
                ../src/vecmat/rmmat.c:  r1matops->ptr_matvecvec=8;
                ../src/vecmat/dlpack.c:  sops->matvecvec=DvechmatVecVec;
                ../src/vecmat/dlpack.c:  sops->ptr_matvecvec=9;
                */
                
                  //printf("File %s line %d A.dsdpops->matvecvec point to %d located in ",__FILE__, __LINE__,A.dsdpops->matvecvec);
                  info=SDPConeVecGetSize(W,&n); DSDPCHKERR(info);
                  info=SDPConeVecGetArray(W,&x); DSDPCHKERR(info);
                  //info=(A.dsdpops->matvecvec)(A.matdata,x,n,v); // DSDPChkDataError(A,info);
                  if(A.dsdpops->ptr_matvecvec == 2){
                    //info=(A.dsdpops->matvecvec)(A.matdata,x,n,v); // DSDPChkDataError(A,info);
                  //static int VechMatVecVec(void* AA, double x[], int n, double *v){
                    {
                      //printf("File %s line %d VechMatVecVec with address %d\n",__FILE__, __LINE__,&VechMatVecVec);
                      vechmat* B=(vechmat*)A.matdata;
                      int info,rank=n,i=0,j,k,kk;
                      const int *ind=B->ind,ishift=B->ishift;
                      double vv=0,dd;
                      const double *val=B->val,nnz=B->nnzeros;
                      
                      if (B->factored==3){
                        info=VechMatGetRank(A.matdata,&rank,n);
                        if (nnz>3 && rank<nnz){
                          info=EigMatVecVec(B->Eig,x,n,&vv);
                          *v=vv*B->alpha;
                          return 0;
                        }
                      }
                      
                      for (k=0; k<nnz; ++k,++ind,++val){
                        kk=*ind-ishift;
                        i=GETI_VECH(kk);
                        j=GETJ_VECH(kk,i);
                        dd=x[i]*x[j]*(*val);
                        vv+=2*dd;
                        if (i==j){ vv-=dd; }
                      }
                      *v=vv*B->alpha;

                      //return 0;
                    }                  
                  
                  }else if(A.dsdpops->ptr_matvecvec == 6){
                  
                    //info=(A.dsdpops->matvecvec)(A.matdata,x,n,v); // DSDPChkDataError(A,info);
                    //static int IdentityMatVecVec(void* AA, double x[], int n, double *v){
                    {  
                      //printf("File %s line %d IdentityMatVecVec with address %d\n",__FILE__, __LINE__,&IdentityMatVecVec);
                      identitymat* B=(identitymat*)A.matdata;
                      int i;
                      *v=0;
                      for (i=0;i<n;i++){
                        *v+=x[i]*x[i];
                      }
                      *v *= B->dm;
                      //return 0;
                    }
                    
                  }else{
                    printf("A.dsdpops->ptr_matvecvec Error!\n");
                  }
                  
                  
                  info=SDPConeVecRestoreArray(W,&x); DSDPCHKERR(info);
                } else {
                  //DSDPNoOperationError(A);
		  exit (-1);
                }
              //}

	      } // end of DSDPDataMatVecVec
              info=DSDPVecAddElement(VAV,vari,aa*aalpha*sum*scl);DSDPCHKVARERR(vari,info);
            }
 	    } // end of else
            //DSDPEventLogEnd(sdpvecvecevent);
          //}  //  original end of DSDPBlockvAv

	  } // end of DSDPBlockvAv
  	} /* End row computations for rank kk of block kk */
   
        }   /* End row computations for all of block kk     */
  
        if (method1==DSDP_TRUE){
  	//info=DSDPBlockADot(&blk[kk].ADATA,1.0,Select,T,MRowI);DSDPCHKBLOCKERR(kk,info);
        //Wei: inline function body of DSDPBlockADot 
        //int DSDPBlockADot(DSDPBlockData *ADATA, double aa, DSDPVec Alpha, DSDPVMat X, DSDPVec AX)
          {
          int    ii,vari,n,nn,info;
          double *x,sum=0,aalpha=0,scl=blk[kk].ADATA.scl;
          //DSDPEventLogBegin(sdpdotevent);   //ignore event
          //info=DSDPVMatScaleDiagonal(T,0.5); //DSDPCHKERR(info);
	  {
          //int DSDPVMatScaleDiagonal(DSDPVMat X, double dscale){
	    DSDPVMat X=T;
	    double dscale=0.5;
            int info;
            if (X.dsdpops->matscalediagonal){
              //printf("X.dsdpops->matscalediagonal %d\n",X.dsdpops->ptr_matscalediagonal);
              // value 1 and 2
              //info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); // DSDPChkMatError(X,info);
              
              /*
                ../src/vecmat/dufull.c:  densematops->matscalediagonal=DTRUMatScaleDiagonal;
                ../src/vecmat/dufull.c:  densematops->ptr_matscalediagonal=1;
                ../src/vecmat/dlpack.c:  densematops->matscalediagonal=DTPUMatScaleDiagonal;
                ../src/vecmat/dlpack.c:  densematops->ptr_matscalediagonal=2;
                */
                
              if(X.dsdpops->ptr_matscalediagonal == 1){
              //info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); // DSDPChkMatError(X,info);
              //static int DTRUMatScaleDiagonal(void* AA, double dd){
              {  
                  dtrumat* B=(dtrumat*) X.matdata;
                  ffinteger LDA=B->LDA;
                  int i,n=B->n;
                  double *v=B->val;
                  for (i=0; i<n; i++){
                    *v*=dscale;
                    v+=LDA+1;    
                  }
                  //return 0;
                }
              }else if(X.dsdpops->ptr_matscalediagonal == 2){
              //info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); // DSDPChkMatError(X,info);
              //static int DTPUMatScaleDiagonal(void* AA, double dd){
              {
                  dtpumat* B=(dtpumat*) X.matdata;
                  int i,n=B->n;
                  double *v=B->val;
                  for (i=0; i<n; i++){
                    *v*=dscale;
                    v+=i+2;    
                  }
                  //return 0;
                }
              
              
              
              }else{
                printf("X.dsdpops->ptr_matscalediagonal Error\n");
              }
              
              
              
            } else {
	      exit (-1);
            }
          //}
	  } // end of DSDPVMatScaleDiagonal 
          //info=DSDPVMatGetSize(T, &n); DSDPCHKERR(info);
	  {
          //int DSDPVMatGetSize(DSDPVMat X,int*n){
	    DSDPVMat X=T;
            int info;
            if (X.dsdpops->matgetsize){
            /*
            ../src/vecmat/dufull.c:  densematops->matgetsize=DTRUMatGetSize;
            ../src/vecmat/dufull.c:  densematops->ptr_matgetsize=1;
            ../src/vecmat/dlpack.c:  densematops->matgetsize=DTPUMatGetSize;
            ../src/vecmat/dlpack.c:  densematops->ptr_matgetsize=2;
            */  
            //printf("X.dsdpops->matgetsize %d\n",X.dsdpops->ptr_matgetsize);
            // 1 and 2
            //info=(X.dsdpops->matgetsize)(X.matdata,&n); //DSDPChkMatError(X,info);
            if(X.dsdpops->ptr_matgetsize == 1){
            //info=(X.dsdpops->matgetsize)(X.matdata,&n); //DSDPChkMatError(X,info);
            //static int DTRUMatGetSize(void *AA, int *n){
            {  
              dtrumat* B=(dtrumat*) X.matdata;
              n=B->n;
              //return 0;
            }
            
            }else if(X.dsdpops->ptr_matgetsize == 2){
            info=(X.dsdpops->matgetsize)(X.matdata,&n); //DSDPChkMatError(X,info);
            //
            //  Need to Fix this, inline cause wrong result!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //
            //
            
            //static int DTPUMatGetSize(void *AA, int *n){
            //{  
              //dtpumat* B=(dtpumat*) X.matdata;
              //n=B->n;
              //return 0;
            //}

            }else{
                printf("X.dsdpops->ptr_matgetsize Error!\n");
            }
            
            
            } else {
              /*
              DSDPNoOperationError(X);
              */
	      exit(-1);
            }
          //}

	  } // end of DSDPVMatGetSize 
          //info=DSDPVMatGetArray(T, &x, &nn); DSDPCHKERR(info);
	  {
          //int DSDPVMatGetArray(DSDPVMat X, double **v, int *nn){
	    DSDPVMat X=T;
	    double **v = &x;
            int info;
            DSDPFunctionBegin;
            if (X.dsdpops->matgeturarray){
            //printf("X.dsdpops->matgeturarray %d\n",X.dsdpops->ptr_matgeturarray);
            /*
            ../src/vecmat/dufull.c:  densematops->matgeturarray=DTRUMatGetDenseArray;
            ../src/vecmat/dufull.c:  densematops->ptr_matgeturarray=1;
            ../src/vecmat/dlpack.c:  densematops->matgeturarray=DTPUMatGetDenseArray;
            ../src/vecmat/dlpack.c:  densematops->ptr_matgeturarray=2;
            */  
            if(X.dsdpops->ptr_matgeturarray == 1){
            //info=(X.dsdpops->matgeturarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
            //static int DTRUMatGetDenseArray(void* A, double *v[], int*n){
                {  
                  dtrumat*  ABA=(dtrumat*)X.matdata;
                  *v=ABA->val;
                  nn=ABA->n*ABA->LDA;
                  //return 0;
                }
                      
            
            }else if(X.dsdpops->ptr_matgeturarray == 2){
                //info=(X.dsdpops->matgeturarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
                //static int DTPUMatGetDenseArray(void* A, double *v[], int*n){
                {  
                  dtpumat*  ABA=(dtpumat*)X.matdata;
                  *v=ABA->val;
                  nn=(ABA->n)*(ABA->n+1)/2;
                  //return 0;
                }
            
            }else{
                printf("X.dsdpops->ptr_matgeturarray Error!\n");
            }
              
            } else {
              *v=0;
              nn=0;
            }
          //}

	  } // end of DSFPVMatGetArray
          for (ii=0;ii<blk[kk].ADATA.nnzmats; ii++){  /* Matrix Entries */
            vari=blk[kk].ADATA.nzmat[ii];
            info=DSDPVecGetElement(Select,vari,&aalpha);DSDPCHKVARERR(vari,info);
            if (aalpha==0.0) continue;
            //info=DSDPDataMatDot(blk[kk].ADATA.A[ii],x,nn,n,&sum);DSDPCHKVARERR(vari,info);
            //Wei: inline function body of DSDPDataMatDot 
            //int DSDPDataMatDot(DSDPDataMat A, double x[], int nn, int n, double *v)
	    {
                int info;
                if (blk[kk].ADATA.A[ii].dsdpops->matdot){
                  //printf("blk[kk].ADATA.A[ii].dsdpops->matdot %d\n",blk[kk].ADATA.A[ii].dsdpops->ptr_matdot);
                  //onlye 2 6 7 and 11.
                  /*
                    ../src/vecmat/drowcol.c:  rcmatoperator->matdot=RCMatDot;
                    ../src/vecmat/drowcol.c:  rcmatoperator->ptr_matdot=1;
                    ../src/vecmat/vech.c:  sops->matdot=VechMatDot;
                    ../src/vecmat/vech.c:  sops->ptr_matdot=2;
                    ../src/vecmat/dufull.c:  sops->matdot=DvecumatDot;
                    ../src/vecmat/dufull.c:  sops->ptr_matdot=3;
                    ../src/vecmat/onemat.c:  cmatops->matdot=ConstMatDot;
                    ../src/vecmat/onemat.c:  cmatops->ptr_matdot=4;
                    ../src/vecmat/zeromat.c:  sops->matdot=ZDot;
                    ../src/vecmat/zeromat.c:  sops->ptr_matdot=5;
                    ../src/vecmat/identity.c:  spdiagops->matdot=IdentityMatDotP;
                    ../src/vecmat/identity.c:  spdiagops->ptr_matdot=6;
                    ../src/vecmat/identity.c:  spdiagops->matdot=IdentityMatDotF;
                    ../src/vecmat/identity.c:  spdiagops->ptr_matdot=11;
                    ../src/vecmat/vechu.c:  sops->matdot=VechMatDot;
                    ../src/vecmat/vechu.c:  sops->ptr_matdot=7;
                    ../src/vecmat/rmmat.c:  r1matops->matdot=R1MatDotP;
                    ../src/vecmat/rmmat.c:  r1matops->ptr_matdot=8;
                    ../src/vecmat/rmmat.c:  r1matops->matdot=R1MatDotU;
                    ../src/vecmat/rmmat.c:  r1matops->ptr_matdot=9;
                    ../src/vecmat/dlpack.c:  sops->matdot=DvechmatDot;
                    ../src/vecmat/dlpack.c:  sops->ptr_matdot=10;
                    */
                  if(blk[kk].ADATA.A[ii].dsdpops->ptr_matdot == 2){
                  //info=(blk[kk].ADATA.A[ii].dsdpops->matdot)(blk[kk].ADATA.A[ii].matdata,x,nn,n,&sum); 
                      //static int VechMatDot(void* AA, double x[], int nn, int n, double *v){
                    {  //printf("File %s line %d VechMatDot with address %d\n",__FILE__, __LINE__,&VechMatDot);
                      vechmat* B=(vechmat*)blk[kk].ADATA.A[ii].matdata;
                      int k,nnz=B->nnzeros; 
                      const int *ind=B->ind;
                      double vv=0, *xx=x-B->ishift;
                      const double *val=B->val;
                      for (k=0;k<nnz;++k,++ind,++val){
                        vv+=(*val)*(*(xx+(*ind)));
                      }
                      sum=2*vv*B->alpha;
                      //return 0;
                    }
                  
                  }else if(blk[kk].ADATA.A[ii].dsdpops->ptr_matdot == 6){
                  //info=(blk[kk].ADATA.A[ii].dsdpops->matdot)(blk[kk].ADATA.A[ii].matdata,x,nn,n,&sum); 
                  //static int IdentityMatDotP(void* AA, double x[], int nn, int n, double *v){
                    {  //printf("File %s line %d IdentityMatDotP with address %d\n",__FILE__, __LINE__,&IdentityMatDotP);
                      identitymat* B=(identitymat*)blk[kk].ADATA.A[ii].matdata;
                      int i;
                      double *xx=x;
                      sum=0;
                      for (i=0;i<n;i++){
                        sum+=*xx;
                        xx+=i+2;
                      }
                      sum *= 2*B->dm;
                      //return 0;
                    }
                  
                  
                  
                  }else if(blk[kk].ADATA.A[ii].dsdpops->ptr_matdot == 7){
                  
                  //static int VechMatDot(void* AA, double x[], int nn, int n, double *v){
                    {
                      //printf("File %s line %d VechMatDot with address %d\n",__FILE__, __LINE__,&VechMatDot);
                      vechmat* B=(vechmat*)blk[kk].ADATA.A[ii].matdata;
                      int k,nnz=B->nnzeros; 
                      const int *ind=B->ind;
                      double vv=0,*xx=x-B->ishift;
                      const double *val=B->val;
                      for (k=0;k<nnz;++k,++ind,++val){
                        vv+=(*val)*(*(xx+(*ind)));
                      }
                      sum=2*vv*B->alpha;
                      //return 0;
                    }
                  
                  }else if(blk[kk].ADATA.A[ii].dsdpops->ptr_matdot == 11){
                  
                  //static int IdentityMatDotF(void* AA, double x[], int nn, int n, double *v){
                    {  //printf("File %s line %d IdentityMatDotF with address %d\n",__FILE__, __LINE__,&IdentityMatDotF);
                      identitymat* B=(identitymat*)blk[kk].ADATA.A[ii].matdata;
                      int i;
                      double *xx=x;
                      sum=0;
                      for (i=0;i<n;i++){
                        sum+=*xx;
                        xx+=n+1;
                      }
                      sum *= 2*B->dm;
                      //return 0;
                    }
                  
                  }else{
                    printf("blk[kk].ADATA.A[ii].dsdpops->ptr_matdot Error!\n");
                  }
                
                  //printf("File %s line %d blk[kk].ADATA.A[ii].dsdpops->matdot point to %d located in ",__FILE__, __LINE__,blk[kk].ADATA.A[ii].dsdpops->matdot);
                  //info=(blk[kk].ADATA.A[ii].dsdpops->matdot)(blk[kk].ADATA.A[ii].matdata,x,nn,n,&sum); 
		  //DSDPChkDataError(blk[kk].ADATA.A[ii],info);
                } else {
                  //DSDPNoOperationError(A);
		  exit(-1);   // error handling omitted
                }

            } //Wei: end of inline function body of DSDPDataMatDot
            info=DSDPVecAddElement(MRowI,vari,1.0*aalpha*sum*scl);DSDPCHKVARERR(vari,info);
          }
          //info=DSDPVMatRestoreArray(T, &x, &nn); DSDPCHKERR(info);
	  {
          //int DSDPVMatRestoreArray(DSDPVMat X, double **v, int *nn){
	    DSDPVMat X=T;
	    double **v = &x;
            int info;
            if (X.dsdpops->matrestoreurarray){
            /*
            ../src/vecmat/dufull.c:  densematops->matrestoreurarray=DTRUMatRestoreDenseArray;
            ../src/vecmat/dufull.c:  densematops->ptr_matrestoreurarray=1;
            ../src/vecmat/dlpack.c:  densematops->matrestoreurarray=DTPUMatRestoreDenseArray;
            ../src/vecmat/dlpack.c:  densematops-ptr_>matrestoreurarray=2;
            */
              if(X.dsdpops->ptr_matrestoreurarray==1){
              //info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
              //static int DTRUMatRestoreDenseArray(void* A, double *v[], int *n){
                {
                  *v=0;nn=0;
                  //return 0;
                }

              }else if(X.dsdpops->ptr_matrestoreurarray==2){
              
              //static int DTPUMatRestoreDenseArray(void* A, double *v[], int *n){
                {
                  *v=0;nn=0;
                  //return 0;
                }
              }else{
                printf("X.dsdpops->ptr_matrestoreurarray Error!\n");
              }
              //info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
            } else {
              *v=0;
              nn=0;
            }
          //}
	  } // end of DSDPVMatRestoreArray
          //info=DSDPVMatScaleDiagonal(T,2.0); DSDPCHKERR(info);
	  {
	  //int DSDPVMatScaleDiagonal(DSDPVMat X, double dscale){
	    DSDPVMat X=T;
	    double dscale=2.0;

	    int info;
	    
	        if (X.dsdpops->matscalediagonal){
              //printf("X.dsdpops->matscalediagonal %d\n",X.dsdpops->ptr_matscalediagonal);
              // value 1 and 2
              //info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); // DSDPChkMatError(X,info);
              
              /*
                ../src/vecmat/dufull.c:  densematops->matscalediagonal=DTRUMatScaleDiagonal;
                ../src/vecmat/dufull.c:  densematops->ptr_matscalediagonal=1;
                ../src/vecmat/dlpack.c:  densematops->matscalediagonal=DTPUMatScaleDiagonal;
                ../src/vecmat/dlpack.c:  densematops->ptr_matscalediagonal=2;
                */
                
              if(X.dsdpops->ptr_matscalediagonal == 1){
              //info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); // DSDPChkMatError(X,info);
              //static int DTRUMatScaleDiagonal(void* AA, double dd){
              {  
                  dtrumat* B=(dtrumat*) X.matdata;
                  ffinteger LDA=B->LDA;
                  int i,n=B->n;
                  double *v=B->val;
                  for (i=0; i<n; i++){
                    *v*=dscale;
                    v+=LDA+1;    
                  }
                  //return 0;
                }
              }else if(X.dsdpops->ptr_matscalediagonal == 2){
              //info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); // DSDPChkMatError(X,info);
              //static int DTPUMatScaleDiagonal(void* AA, double dd){
              {
                  dtpumat* B=(dtpumat*) X.matdata;
                  int i,n=B->n;
                  double *v=B->val;
                  for (i=0; i<n; i++){
                    *v*=dscale;
                    v+=i+2;    
                  }
                  //return 0;
                }
              
              
              
              }else{
                printf("X.dsdpops->ptr_matscalediagonal Error\n");
              }
              
              
              
            }
            else {
	      exit (-1);
	    }
	  //}
	  }  // end of DSDPVMatScaleDiagonal

          //DSDPEventLogEnd(sdpdotevent);
          }   //Wei: end of inline function body of DSDPBlockADot
        }   /* End row computations for all of block ll     */
      }     /* End row computations for all blocks          */
      info=DSDPVecAddElement(vrhs1,i,rhs1i);DSDPCHKERR(info);
      info=DSDPVecAddElement(vrhs2,i,rhs2i);DSDPCHKERR(info);
      //info=DSDPSchurMatAddRow(M,i,1.0,MRowI);DSDPCHKERR(info);
      {
      //int DSDPSchurMatAddRow(DSDPSchurMat M, int row, double alpha, DSDPVec R){
	int row = i;
	double alpha=1.0;
	DSDPVec R = MRowI;
        int info,j,m;
        double *v,rr,dd=1e-1*M.schur->dd;
        DSDPVec rhs3=M.schur->rhs3;
        DSDPTruth flag;
        info=DSDPVecGetSize(R,&m); DSDPCHKERR(info);
        if (row==0){
        } else if (row==m-1){
          info=DSDPVecGetR(R,&rr);DSDPCHKERR(info);
          info=DSDPVecAddR(rhs3,alpha*rr);DSDPCHKERR(info);
        } else if (M.dsdpops->mataddrow){
          info=DSDPVecGetArray(R,&v); DSDPCHKERR(info);
          /*    v[row]=DSDPMax(0,v[row]); v[row]+=1.0e-15; */
          for (j=0;j<m;j++){ if (fabs(v[j]) < 1e-25 && row!=j){ v[j]=0.0;} }
          v[row]*=(1.0+dd);
          //info=DSDPZeroFixedVariables(M,R);DSDPCHKERR(info);
	  {
          //int DSDPZeroFixedVariables( DSDPSchurMat M, DSDPVec dy){
            int i,info; 
            FixedVariables *fv=&M.schur->fv;
            for (i=0;i<fv->nvars;i++){
              info=DSDPVecSetElement(R,fv->var[i],0.0);DSDPCHKERR(info);
            }
          //} end of original end

          } // end of DSDPZeroFixedVariables

          //info=DSDPIsFixed(M,row,&flag);DSDPCHKERR(info); 
	  {
	  //int DSDPIsFixed( DSDPSchurMat M, int vari, DSDPTruth *flag){
	    int vari = row;
	    int i;
	    FixedVariables *fv=&M.schur->fv;
	    flag=DSDP_FALSE;
	    for (i=0;i<fv->nvars;i++){
	      if (fv->var[i]==vari){
	        flag=DSDP_TRUE;
	        break;
	      }
	    }
	  //} original end 

	  } // end of DSDPIsFixed 

          if (flag==DSDP_TRUE){info=DSDPVecSetBasis(R,row);DSDPCHKERR(info);}
          
          //printf("M.dsdpops->mataddrow %d\n",M.dsdpops->ptr_mataddrow);
          // only 1
          /*
            ../src/vecmat/dufull.c:  mops->mataddrow=DTRUMatAddRow;
            ../src/vecmat/dufull.c:  mops->ptr_mataddrow=1;
            ../src/vecmat/cholmat.c:  mops->mataddrow=Taddline;
            ../src/vecmat/cholmat.c:  mops->ptr_mataddrow=2;
            ../src/vecmat/diag.c:  sops->mataddrow=DiagMatAddRow2;
            ../src/vecmat/diag.c:  sops->ptr_mataddrow=3;
            ../src/vecmat/dlpack.c:  mops->mataddrow=DTPUMatAddRow;
            ../src/vecmat/dlpack.c:  mops->ptr_mataddrow=4;
            */
          if(M.dsdpops->ptr_mataddrow == 1){
            //info=(M.dsdpops->mataddrow)(M.data,row-1,alpha,v+1,m-2); // DSDPChkMatError(M,info);
            //static int DTRUMatAddRow(void* AA, int nrow, double dd, double row[], int n){
            {  
              dtrumat* B=(dtrumat*) M.data;
              ffinteger ione=1,LDA=B->LDA,nn,INCX=1,INCY=B->LDA;
              double *vv=B->val;
              int nrow=row-1;
              nn=nrow;
              daxpy(&nn,&alpha,v+1,&INCX,vv+nrow,&INCY);
              nn=nrow+1;
              daxpy(&nn,&alpha,v+1,&ione,vv+nrow*LDA,&ione);

              //return 0;
            }
          }else{
            printf("M.dsdpops->ptr_mataddrow Error!\n");
          }
          
          //info=(M.dsdpops->mataddrow)(M.data,row-1,alpha,v+1,m-2); // DSDPChkMatError(M,info);
          info=DSDPVecRestoreArray(R,&v); DSDPCHKERR(info);  
          info=DSDPVecGetR(R,&rr); DSDPCHKERR(info);
          info=DSDPVecAddElement(rhs3,row,alpha*rr); DSDPCHKERR(info);
        } else {
          //DSDPNoOperationError(M);
	  exit (-1);
        }
      //} end of original DSDPSchurMatAddRow
      } // end of DSDPSchurMatAddRow
    }
  //} // end of original SDPConeComputeHessian
  }   //Wei: end of inlined body of SDPConeComputeHessian
  DSDPFunctionReturn(0);   
}

#undef __FUNCT__  
#define __FUNCT__ "KDPConeMultiply"
static int KSDPConeMultiply( void *K, double mu, DSDPVec vrow, DSDPVec vin, DSDPVec vout){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  SDPConeValid(sdpcone);
  DSDPFunctionBegin;
  for (kk=0; kk<sdpcone->nblocks; kk++){      
    info=SDPConeMultiply(sdpcone,kk,mu,vrow,vin,vout);DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0);   
}

#undef __FUNCT__  
#define __FUNCT__ "KDPConeRHS  "
static int KSDPConeRHS( void *K, double mu, DSDPVec vrow, DSDPVec vrhs1, DSDPVec vrhs2){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0; kk<sdpcone->nblocks; kk++){      
    if (sdpcone->blk[kk].n<1) continue;
    //info=SDPConeComputeRHS(sdpcone,kk,mu,vrow,vrhs1,vrhs2); DSDPCHKBLOCKERR(kk,info);
    //Wei: inlined function body of SDPConeComputeRHS
    {
    //int SDPConeComputeRHS( SDPCone sdpcone, int blockj, double mu, DSDPVec vrow, DSDPVec vrhs1, DSDPVec vrhs2){
      int blockj = kk; 
      int info,i,ii,k,rank,nnzmats;
      double dtmp,dyiscale=1,ack,scl,rtemp;
      SDPblk *sdp=&sdpcone->blk[blockj];
      SDPConeVec W=sdp->W,W2=sdp->W2;
      DSDPDataMat AA;
      DSDPVMat T=sdp->T;
      DSDPDualMat S=sdp->S;
      DSDPTruth method1;
      DSDPIndex IS=sdp->IS;
    
      info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
      method1=DSDP_TRUE;
      method1=DSDP_FALSE;
    
      if (method1==DSDP_TRUE){
        info=DSDPBlockCountNonzeroMatrices(&sdp->ADATA,&nnzmats);DSDPCHKERR(info);
        for (i=0;i<nnzmats;i++){
          info=DSDPBlockGetMatrix(&sdp->ADATA,i,&ii,&scl,&AA);DSDPCHKERR(info);
          info=DSDPVecGetElement(vrow,ii,&dyiscale);DSDPCHKVARERR(ii,info);
          if (dyiscale==0) continue;
          info=DSDPDataMatGetRank(AA,&rank,sdp->n);DSDPCHKVARERR(ii,info);
          for (k=0; k<rank; k++){
    	info=DSDPDataMatGetEig(AA,k,W,IS,&ack); DSDPCHKVARERR(ii,info);
    	if (ack==0) continue;
    	//info=DSDPDualMatInverseMultiply(S,IS,W,W2);DSDPCHKVARERR(ii,info);
	{
        //int DSDPDualMatInverseMultiply(DSDPDualMat S, DSDPIndex IS, SDPConeVec B, SDPConeVec X){
	  SDPConeVec B=W;
          SDPConeVec X=W2;
          int info,n;
          double *bb,*xx;
          //DSDPEventLogBegin(sdpdualsolve);
          printf("Executing S.dsdpops->matinversemultiply\n");
          //
          //    This section does not execute...................double check..................!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //
          //
          if (S.dsdpops->matinversemultiply){
            printf("S.dsdpops->matinversemultiply %d\n",S.dsdpops->ptr_matinversemultiply);
          
          
          
          
            info=SDPConeVecGetSize(X,&n); DSDPCHKERR(info);
            info=SDPConeVecGetArray(B,&bb); DSDPCHKERR(info);
            info=SDPConeVecGetArray(X,&xx); DSDPCHKERR(info);
            info=(S.dsdpops->matinversemultiply)(S.matdata,IS.indx+1,IS.indx[0],bb,xx,n); DSDPChkDMatError(S,info);
            info=SDPConeVecRestoreArray(X,&xx); DSDPCHKERR(info);
            info=SDPConeVecRestoreArray(B,&bb); DSDPCHKERR(info);
          } else {
            //DSDPNoOperationError(S);
	    exit (-1);
          }
          //DSDPEventLogEnd(sdpdualsolve);
        //}  end of original

	} // end of DSDPDualMatInverseMultiply
    	//info=SDPConeVecDot(W,W2,&rtemp); DSDPCHKVARERR(ii,info);
	{
        //int SDPConeVecDot(SDPConeVec V1, SDPConeVec V2, double *ans){
	  double *ans = &rtemp;
	  SDPConeVec V1 = W;
	  SDPConeVec V2 = W2;
          ffinteger ione=1, nn=V1.dim;
          double *v1=V1.val,*v2=V2.val;
          *ans=ddot(&nn,v1,&ione,v2,&ione);
        //}
	} // end of SDPConeVecDot
    	dtmp=rtemp*ack*mu*dyiscale*scl;
    	info=DSDPVecAddElement(vrhs2,ii,dtmp);DSDPCHKVARERR(ii,info);
          }
        }
        
      } else {
        //info=DSDPVMatZeroEntries(T);DSDPCHKERR(info);
        {
        //int DSDPVMatZeroEntries(DSDPVMat X){
	  DSDPVMat X=T; 
          int info;
          
          if (X.dsdpops->matzeroentries){
            //info=(X.dsdpops->matzeroentries)(X.matdata); //DSDPChkMatError(X,info);
            //printf("X.dsdpops->matzeroentries %d\n",X.dsdpops->ptr_matzeroentries);
            // only 1 and 4
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
            if(X.dsdpops->ptr_matzeroentries==1){
            //info=(X.dsdpops->matzeroentries)(X.matdata); //DSDPChkMatError(X,info);
            //static int DTRUMatZero(void* AA){
            {
              //printf("File %s line %d DTRUMatZero with address %d\n",__FILE__, __LINE__,&DTRUMatZero);
              dtrumat* B=(dtrumat*) X.matdata;
              int mn=B->n*(B->LDA);
              double *vv=B->val;
              memset((void*)vv,0,mn*sizeof(double));
              B->status=Assemble;
              //return 0;
            }
            }else if(X.dsdpops->ptr_matzeroentries==4){
              //printf("File %s line %d X.dsdpops->matzeroentries point to %d located in ",__FILE__, __LINE__,X.dsdpops->matzeroentries);
              //info=(X.dsdpops->matzeroentries)(X.matdata); //DSDPChkMatError(X,info);
            //static int DTPUMatZero(void* AA){
            {
              //printf("File %s line %d DTPUMatZero with address %d\n",__FILE__, __LINE__,&DTPUMatZero);
              dtpumat* B=(dtpumat*) X.matdata;
              int mn=B->n*(B->n+1)/2;
              double *vv=B->val;
              memset((void*)vv,0,mn*sizeof(double));
              //return 0;
            }
            }else{
                printf("X.dsdpops->ptr_matzeroentries Error! Should be impossible to get here!\n");
            }
              
              
            }else {
	    exit(-1);
            //DSDPNoOperationError(X);
          }
        //}
	} // end of inlining DSDPVMatZeroEntries
        info=DSDPDualMatInverseAdd(S,mu,T);DSDPCHKERR(info);
        //info=DSDPBlockADot(&sdp->ADATA,1.0,vrow,T,vrhs2);DSDPCHKERR(info);
        //Wei: inlined function body of DSDPBlockADot
        {
        //int DSDPBlockADot(DSDPBlockData *ADATA, double aa, DSDPVec Alpha, DSDPVMat X, DSDPVec AX){
          DSDPBlockData *ADATA = &sdp->ADATA;
          double aa = 1.0;
          DSDPVec Alpha = vrow;
          DSDPVMat X=T; 
          DSDPVec AX=vrhs2;
    
          int    ii,vari,n,nn,info;
          double *x,sum=0,aalpha=0,scl=ADATA->scl;
        
          //DSDPEventLogBegin(sdpdotevent);
          info=DSDPVMatScaleDiagonal(X,0.5); DSDPCHKERR(info);
          info=DSDPVMatGetSize(X, &n); DSDPCHKERR(info);
          info=DSDPVMatGetArray(X, &x, &nn); DSDPCHKERR(info);
          for (ii=0;ii<ADATA->nnzmats; ii++){  /* Matrix Entries */
            vari=ADATA->nzmat[ii];
            info=DSDPVecGetElement(Alpha,vari,&aalpha);DSDPCHKVARERR(vari,info);
            if (aalpha==0.0) continue;
            info=DSDPDataMatDot(ADATA->A[ii],x,nn,n,&sum);DSDPCHKVARERR(vari,info);
            info=DSDPVecAddElement(AX,vari,aa*aalpha*sum*scl);DSDPCHKVARERR(vari,info);
          }
          info=DSDPVMatRestoreArray(X, &x, &nn); DSDPCHKERR(info);
	  {
          //int DSDPVMatRestoreArray(DSDPVMat X, double **v, int *nn){
	    double **v = &x;
            int info;
            
            if (X.dsdpops->matrestoreurarray){
            /*
            ../src/vecmat/dufull.c:  densematops->matrestoreurarray=DTRUMatRestoreDenseArray;
            ../src/vecmat/dufull.c:  densematops->ptr_matrestoreurarray=1;
            ../src/vecmat/dlpack.c:  densematops->matrestoreurarray=DTPUMatRestoreDenseArray;
            ../src/vecmat/dlpack.c:  densematops-ptr_>matrestoreurarray=2;
            */
              if(X.dsdpops->ptr_matrestoreurarray==1){
              //info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
              //static int DTRUMatRestoreDenseArray(void* A, double *v[], int *n){
                {
                  *v=0;nn=0;
                  //return 0;
                }

              }else if(X.dsdpops->ptr_matrestoreurarray==2){
              
              //static int DTPUMatRestoreDenseArray(void* A, double *v[], int *n){
                {
                  *v=0;nn=0;
                  //return 0;
                }
              }else{
                printf("X.dsdpops->ptr_matrestoreurarray Error!\n");
              }
              //info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
            } else {
              *v=0;
              nn=0;
            }
          //}
	  } // end of DSDPVMatRestoreArray
          info=DSDPVMatScaleDiagonal(X,2.0); DSDPCHKERR(info);
          //DSDPEventLogEnd(sdpdotevent);
        //}
    
        } // end of inlined function body of DSDPBlockADot
      }
    //}
    } // end of inlined function body of SDPConeComputeRHS
  } 
  DSDPFunctionReturn(0);   
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSetup"
static int KSDPConeSetup(void* K, DSDPVec y){
  int info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=SDPConeSetup(sdpcone,y);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSetup2"
static int KSDPConeSetup2(void* K, DSDPVec y, DSDPSchurMat M){
  int info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  info=SDPConeSetup2(sdpcone,y,M); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeDestroy"
static int KSDPConeDestroy(void* K){
  int info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=SDPConeDestroy(sdpcone);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSize"
static int KSDPConeSize(void* K,double *n){
  SDPCone sdpcone=(SDPCone)K;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  *n=sdpcone->nn;
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSparsity"
static int KSDPConeSparsity(void *K,int row, int *tnnz, int rnnz[], int m){
  SDPCone sdpcone=(SDPCone)K;
  SDPblk *blk=sdpcone->blk;
  int info,j,kk;
  int nnzblocks=sdpcone->ATR.nnzblocks[row],*nzblocks=sdpcone->ATR.nzblocks[row];
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (j=0; j<nnzblocks; j++){
    kk=nzblocks[j];
    if (blk[kk].n<1) continue;
    info=DSDPBlockDataMarkNonzeroMatrices(&blk[kk].ADATA,rnnz);DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "KSDPConeComputeSS"
static int KSDPConeComputeSS(void *K, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  SDPblk *blk=sdpcone->blk;
  DSDPTruth psdefinite;
  DSDPDualMat SS;

  DSDPFunctionBegin;
  psdefinite = DSDP_TRUE;
  for (kk=sdpcone->nblocks-1; kk>=0 && psdefinite == DSDP_TRUE; kk--){
    if (blk[kk].n<1) continue;
    if (flag==DUAL_FACTOR) SS=blk[kk].S;
    else SS=blk[kk].SS;

    switch (sdpcone->optype){
    default:
      //info=SDPConeComputeSS(sdpcone,kk,Y,blk[kk].T);DSDPCHKBLOCKERR(kk,info);
      //We: inline function body of SDPConeComputeSS
      { 
      //int SDPConeComputeSS(SDPCone sdpcone, int blockj, DSDPVec Y, DSDPVMat SS){
	int blockj = kk;
        DSDPVMat SS=blk[kk].T;

        int info;
        DSDPFunctionBegin;
        info=DSDPVMatZeroEntries(SS); DSDPCHKBLOCKERR(blockj,info);
        //info=DSDPBlockASum(&sdpcone->blk[blockj].ADATA,1,Y,SS); DSDPCHKBLOCKERR(blockj,info);
        //Wei: inlined function body of DSDPBlockASum
        {
        //int DSDPBlockASum(DSDPBlockData *ADATA, double aa, DSDPVec Yk, DSDPVMat XX){
          DSDPBlockData *ADATA=&sdpcone->blk[blockj];
          double aa = 1;
          DSDPVec Yk = Y;
          DSDPVMat XX = SS;
       
          double *xx,ytmp,scl=ADATA->scl;
          int    ii,vari,n,nn,info;
        
          info=DSDPVMatGetSize(XX, &n); DSDPCHKERR(info);
          info=DSDPVMatGetArray(XX, &xx, &nn); DSDPCHKERR(info);
          for (ii=0;ii<ADATA->nnzmats;ii++){
            vari=ADATA->nzmat[ii];
            info=DSDPVecGetElement(Yk,vari,&ytmp);DSDPCHKVARERR(vari,info);
            if (ytmp==0) continue;
            info = DSDPDataMatAddMultiple(ADATA->A[ii], -aa*scl*ytmp, xx,nn,n); DSDPCHKVARERR(vari,info);
          }
          //info=DSDPVMatRestoreArray(XX, &xx, &nn); DSDPCHKERR(info);
	  {
	  //int DSDPVMatRestoreArray(DSDPVMat X, double **v, int *nn){
	    DSDPVMat X=XX;
	    double **v=&xx;
	    
	    int info;
	    //DSDPFunctionBegin;
	     
	    if (X.dsdpops->matrestoreurarray){
            /*
            ../src/vecmat/dufull.c:  densematops->matrestoreurarray=DTRUMatRestoreDenseArray;
            ../src/vecmat/dufull.c:  densematops->ptr_matrestoreurarray=1;
            ../src/vecmat/dlpack.c:  densematops->matrestoreurarray=DTPUMatRestoreDenseArray;
            ../src/vecmat/dlpack.c:  densematops-ptr_>matrestoreurarray=2;
            */
              if(X.dsdpops->ptr_matrestoreurarray==1){
              //info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
              //static int DTRUMatRestoreDenseArray(void* A, double *v[], int *n){
                {
                  *v=0;nn=0;
                  //return 0;
                }

              }else if(X.dsdpops->ptr_matrestoreurarray==2){
              
              //static int DTPUMatRestoreDenseArray(void* A, double *v[], int *n){
                {
                  *v=0;nn=0;
                  //return 0;
                }
              }else{
                printf("X.dsdpops->ptr_matrestoreurarray Error!\n");
              }
              //info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
            }else {
	      *v=0;
	      nn=0;
	    }
	    //DSDPFunctionReturn(0);
	  //} end of original function body
	  } // end of DSDPVMatRestoreArray
        //} //original end before inlining
        }  // end of inlined function body of DSDPBlockASum
      //}  // end of original SDPConeCOmputeSS before inlining
      } // end of SDPConeComputeSS
      info=DSDPDualMatSetArray(SS,blk[kk].T); DSDPCHKBLOCKERR(kk,info);
      info=DSDPDualMatCholeskyFactor(SS,&psdefinite); DSDPCHKBLOCKERR(kk,info);
      if (psdefinite == DSDP_FALSE){
	if (flag==DUAL_FACTOR){
	  DSDPLogInfo(0,2,"Dual SDP Block %2.0f not PSD\n",kk);
	} else {
	  DSDPLogInfo(0,2,"Primal SDP Block %2.0f not PSD\n",kk);
	}
      }
      break;
    }
  }
  *ispsdefinite=psdefinite;
  if (flag==DUAL_FACTOR && psdefinite==DSDP_TRUE){
    info=DSDPVecCopy(Y,sdpcone->YY);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSDPConeInvertSS"
static int KSDPConeInvertSS(void *K){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  DSDPDualMat SS;
  
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0;kk<sdpcone->nblocks;kk++){
    if (sdpcone->blk[kk].n<1) continue;
    SS=sdpcone->blk[kk].S;
    info=DSDPDualMatInvert(SS);DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "KSDPConeComputeMaxStepLength"
static int KSDPConeComputeMaxStepLength(void *K, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
  int kk,info;
  double smaxstep,maxmaxstep=1.0e20;
  SDPCone sdpcone=(SDPCone)K;
  DSDPDualMat SS;
  SDPblk *blk=sdpcone->blk;
  DSDPDSMat DS;
  DSDPVMat T;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0; kk<sdpcone->nblocks; kk++){      
    if (blk[kk].n<1) continue;
    if (flag==DUAL_FACTOR) SS=blk[kk].S;
    else SS=blk[kk].SS;
    DS=blk[kk].DS; T=blk[kk].T;

    info=SDPConeComputeSS(sdpcone,kk,DY,T);DSDPCHKBLOCKERR(kk,info);
    info=DSDPDSMatSetArray(DS,T); DSDPCHKBLOCKERR(kk,info);

    //info=DSDPLanczosStepSize( &blk[kk].Lanczos,blk[kk].W,blk[kk].W2,SS,DS,&smaxstep );DSDPCHKBLOCKERR(kk,info);
    //Wei: inlined function body of DSDPLanczosStepSize
    {
    //int DSDPLanczosStepSize( DSDPLanczosStepLength *LZ, SDPConeVec W1, SDPConeVec W2, DSDPDualMat S, DSDPDSMat DS, double *maxstep ){
      DSDPLanczosStepLength *LZ = &blk[kk].Lanczos ;    // parameter 1 passing
      SDPConeVec W1 = blk[kk].W;                        // parameter 2 passing
      SDPConeVec W2 = blk[kk].W2;                       // parameter 3 passing
      DSDPDualMat S = SS;                               // parameter 4 passing, 5 the same
      double*  maxstep = &smaxstep;                     // parameter 6 passing

      
      int info,m;
      double smaxstep,mineig;
      struct _P_Mat3 PP;
      Mat3 A=&PP; 
    
      A->ss=S;
      A->ds=DS; A->V=W2;
      A->type=1;
      m=LZ->lanczosm;
    
      if (LZ->type==1){
        //info = ComputeStepFAST(A,LZ->Q,m,W1,LZ->dwork4n,LZ->iwork10n,&smaxstep,&mineig);DSDPCHKERR(info);
        //Wei: inlined function body of ComputeStepFAST
        {
        //static int ComputeStepFAST(Mat3 A, SDPConeVec *Q, int m, SDPConeVec W, double *dwork, int *iwork,double *maxstep ,double *mineig){
        
          int i,j,n,info;
          double tt,wnorm, phi;
          double lambda1=0,lambda2=0,delta=0;
          double res1,res2,beta;
          double one=1.0;
          int N=m;
          double *diag,*subdiag,*ddwork;
        
          diag=LZ->dwork4n;
          subdiag=LZ->dwork4n+m;
          ddwork=LZ->dwork4n+2*m;
        
          if (A->type==1){
            for (i=0; i<m; i++){ diag[i]=-1; subdiag[i]=0;}
          } else {
            for (i=0; i<m; i++){ diag[i]=1.0; subdiag[i]=0;}
          }
          info = SDPConeVecSet(one,LZ->Q[0]);DSDPCHKERR(info);
          info = SDPConeVecNormalize(LZ->Q[0]);DSDPCHKERR(info);
        
          for (i=0; i<m; i++){
            //info = MatMult3(A,LZ->Q[0],W1);DSDPCHKERR(info);
	    //Wei: inline function body of MatMult3
            {
	    //int MatMult3(Mat3 A, SDPConeVec X, SDPConeVec Y){
	      SDPConeVec X = LZ->Q[0];
	      SDPConeVec Y = W1; 
	      int info=0;
	      double minus_one=-1.0;
	    
	      /*  DSDPEventLogBegin(id2); */
	      if (A->type==2){
	        info=DSDPVMatMult(A->x,X,Y);DSDPCHKERR(info);
	        //printf("A->type==2\n");
	      } else {
	        //printf("A->type!=2\n");
	        info=DSDPDualMatCholeskySolveBackward(A->ss,X,Y); DSDPCHKERR(info);
	        info=DSDPDSMatMult(A->ds,Y,A->V); DSDPCHKERR(info);
	        info=DSDPDualMatCholeskySolveForward(A->ss,A->V,Y); DSDPCHKERR(info);
	        info=SDPConeVecScale(minus_one,Y); DSDPCHKERR(info);
	      }
	      /*  DSDPEventLogEnd(id2);*/
	    //} // original end of MatMult3
	    
            } //end of inline function body MatMult3
            info = SDPConeVecNorm2(W1,&phi);DSDPCHKERR(info);
            if (phi!=phi){ smaxstep = 0.0;  return 0;} 
            if (i>0){
              tt=-subdiag[i-1];
              info = SDPConeVecAXPY(tt,LZ->Q[1],W1);DSDPCHKERR(info);
            }
            info = SDPConeVecDot(W1,LZ->Q[0],&tt);DSDPCHKERR(info);
            diag[i]=tt;
            tt*=-1.0;
            info = SDPConeVecAXPY(tt,LZ->Q[0],W1);DSDPCHKERR(info);
            info = SDPConeVecNorm2(W1,&wnorm);DSDPCHKERR(info);
            if (wnorm <= 1.0 * phi){
              for (j=0;j<=i;j++){
        	if (j==i-1){
        	  info = SDPConeVecDot(W1,LZ->Q[1],&tt);DSDPCHKERR(info);
        	  if (tt==tt){tt*=-1.0;} else {tt=0;}
        	  info = SDPConeVecAXPY(tt,LZ->Q[1],W1);DSDPCHKERR(info);
        	  subdiag[i-1]-=tt;
        	} else if (j==i){
        	  info = SDPConeVecDot(W1,LZ->Q[0],&tt);DSDPCHKERR(info);
        	  if (tt==tt){tt*=-1.0;} else {tt=0;}
        	  info = SDPConeVecAXPY(tt,LZ->Q[0],W1);DSDPCHKERR(info);
        	  diag[i]-=tt; 
        	}
        
              } 
            }
            
            info = SDPConeVecNorm2(W1,&wnorm);DSDPCHKERR(info);
            /*    printf("PHI: %4.4e, VNORM: %4.2e Diag: %4.2e\n",phi,wnorm,diag[i]); */
            if (i<m-1){
              subdiag[i]=wnorm;
            }
            if (fabs(wnorm)<=1.0e-10){i++;break;}
            info=SDPConeVecCopy(LZ->Q[0],LZ->Q[1]);DSDPCHKERR(info);
            info=SDPConeVecCopy(W1,LZ->Q[0]);DSDPCHKERR(info);
            info=SDPConeVecNormalize(LZ->Q[0]); DSDPCHKERR(info);
            
          }
          
          /*  DSDPEventLogBegin(id1); */
          info=DSDPGetTriDiagonalEigs(m,diag,subdiag,ddwork,LZ->iwork10n);  DSDPCHKERR(info);
          /*  DSDPEventLogEnd(id1); */
          if (N==0){
            lambda1=-0.0;
            delta=1.0e-20;
            mineig=0;
          } else if (N==1){
            lambda1=-diag[0];
            delta=1.0e-20;
            mineig=diag[0];
          } else if (N>1){
            lambda1=-diag[N-1];
            lambda2=-diag[N-2];
        
            res1=1.0e-8;
            res2=1.0e-8;
             
            tt = -lambda1 + lambda2 - res2;
            if (tt>0) beta=tt;
            else beta=1.0e-20;
            delta = DSDPMin(res1,sqrt(res1)/beta);
            
            mineig=diag[0];
          }
        
          
          if (delta-lambda1>0)
            smaxstep = 1.0/(delta-lambda1);
          else
            smaxstep = 1.0e+30;
        
          info=SDPConeVecGetSize(W1,&n);DSDPCHKERR(info);
          DSDPLogInfo(0,19,"Step Length: Fast Lanczos Iterates: %2d, Max: %d, Block Size: %d, VNorm: %3.1e, Lambda1: %4.4e, Lambda2: %4.4e, Delta: %4.2e, Maxstep: %4.2e\n",
        	      i,m,n,wnorm,lambda1,lambda2,delta,smaxstep);
        
        //}
    
    
    
        } // end of inlined function body of COmputeStepFAST
        *maxstep=smaxstep;
      } else if (LZ->type==2){
        info = ComputeStepROBUST(A,LZ->Q,m,LZ->Q[m],W1,LZ->darray/*LZ->TT*/,LZ->Tv,LZ->dwork4n,&smaxstep,&mineig);DSDPCHKERR(info);
        *maxstep=smaxstep;
      } else {
        DSDPSETERR1(1,"Lanczos Step Length Has not been SetUp. Type: %d\n",LZ->type);
      }
    //}


    } // end of inlined function body of DSDPLanczosStepSize
    DSDPLogInfo(0,12,"Block %d, PD %d, maxstepsize: %4.4e\n",kk,flag,smaxstep);
    maxmaxstep=DSDPMin(smaxstep,maxmaxstep); 
  }
  *maxsteplength=maxmaxstep;
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "KSDPConeAddANorm2"
static int KSDPConeAddANorm2(void *K, DSDPVec ANorm2){
  int kk,info;
  SDPCone sdpcone=(SDPCone)K;
  SDPblk *blk=sdpcone->blk;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0; kk<sdpcone->nblocks; kk++){      
    if (blk[kk].n<1) continue;
    info=DSDPBlockANorm2( &blk[kk].ADATA,ANorm2,blk[kk].n); DSDPCHKBLOCKERR(kk,info);
  }
  DSDPFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "KSDPConeSetX"
static int KSDPConeSetX(void *K, double mu, DSDPVec Y,DSDPVec DY){
  SDPCone sdpcone=(SDPCone)K;
  int info;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=DSDPVecCopy(Y,sdpcone->YX);DSDPCHKERR(info);
  info=DSDPVecCopy(DY,sdpcone->DYX);DSDPCHKERR(info);
  sdpcone->xmakermu=mu;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeComputeXX"
static int KSDPConeComputeXX(void *K, double mu, DSDPVec Y,DSDPVec DY,DSDPVec AX,double* tracexs){

  SDPCone sdpcone=(SDPCone)K;
  int info,kk;
  double xnorm,trxs,xtrace;
  DSDPVMat X;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=KSDPConeSetX(K,mu,Y,DY);DSDPCHKERR(info);
  for (kk=0; kk<sdpcone->nblocks; kk++){
    if (sdpcone->blk[kk].n<1) continue;
    X=sdpcone->blk[kk].T;
    info=SDPConeComputeX3(sdpcone,kk,mu,Y,DY,X);DSDPCHKBLOCKERR(kk,info);
    info=SDPConeComputeXDot(sdpcone,kk,Y,X,AX,&xtrace,&xnorm,&trxs);DSDPCHKBLOCKERR(kk,info);
    *tracexs+=trxs;
    DSDPLogInfo(0,10,"SDP Block %d: norm(X): %4.2e, trace(X): %4.2e, trace(XS): %4.2e.\n",kk,xnorm,xtrace,trxs);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "KSDPConeComputeLogSDeterminant"
static int KSDPConeComputeLogSDeterminant(void *K, double *logdetobj, double *logdet){
  int kk,info;
  double dlogdet=0,dlogdet2=0,dd;
  SDPCone sdpcone=(SDPCone)K;
  SDPblk *blk=sdpcone->blk;

  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  for (kk=0; kk<sdpcone->nblocks; kk++){
    if (blk[kk].n<1) continue;
    info=DSDPDualMatLogDeterminant(blk[kk].S,&dd);DSDPCHKBLOCKERR(kk,info);
    dlogdet+=dd*blk[kk].gammamu;
    dlogdet2+=dd*blk[kk].bmu;
  }
  *logdet=dlogdet;
  *logdetobj=dlogdet2;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSDPConeMonitor"
int KSDPConeMonitor(void *K, int tag){
  DSDPFunctionBegin;
  DSDPFunctionReturn(0);
}

static struct DSDPCone_Ops kops;
static const char *sdpconename ="SDP Cone";

#undef __FUNCT__
#define __FUNCT__ "SDPConeOperationsInitialize"
static int SDPConeOperationsInitialize(struct  DSDPCone_Ops* coneops){
  int info;
  if (coneops==NULL) return 0;
  info=DSDPConeOpsInitialize(coneops); DSDPCHKERR(info);
  coneops->conehessian=KSDPConeComputeHessian;
  coneops->conerhs=KSDPConeRHS;
  coneops->conesetup=KSDPConeSetup;
  coneops->conesetup2=KSDPConeSetup2;
  coneops->conedestroy=KSDPConeDestroy;
  coneops->conecomputes=KSDPConeComputeSS;
  coneops->coneinverts=KSDPConeInvertSS;
  coneops->conesetxmaker=KSDPConeSetX;
  coneops->conecomputex=KSDPConeComputeXX;
  coneops->conemaxsteplength=KSDPConeComputeMaxStepLength;
  coneops->conelogpotential=KSDPConeComputeLogSDeterminant;
  coneops->conesize=KSDPConeSize;
  coneops->conesparsity=KSDPConeSparsity;
  coneops->conehmultiplyadd=KSDPConeMultiply;
  coneops->coneanorm2=KSDPConeAddANorm2;
  coneops->conemonitor=KSDPConeMonitor;
  coneops->id=1;
  coneops->name=sdpconename;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPAddSDP"
/*!
\fn int DSDPAddSDP(DSDP dsdp, SDPCone sdpcone);
\brief Pass a semidefinite cone to the solver.
\param dsdp solver
\param sdpcone semidefinite cone
*/
int DSDPAddSDP(DSDP dsdp,SDPCone sdpcone){
  int info;
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=SDPConeOperationsInitialize(&kops); DSDPCHKERR(info);
  info=DSDPAddCone(dsdp,&kops,(void*)sdpcone); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

