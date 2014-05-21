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
      info=DSDPVecZero(MRowI);DSDPCHKERR(info);
      //info=DSDPSchurMatRowColumnScaling(M,i,Select,&ncols); DSDPCHKERR(info); 
      {
      //int DSDPSchurMatRowColumnScaling(DSDPSchurMat M,int row, DSDPVec V, int *nzcols){
	int row = i;
        DSDPVec V=Select;
        int *nzcols = &ncols;

        int info,m;
        double *cols,r=M.schur->r;
        DSDPTruth flag;
        info=DSDPVecSet(0.0,V);DSDPCHKERR(info);
        info=DSDPVecGetSize(V,&m);DSDPCHKERR(info);
        if (row==0){info=DSDPVecZero(V);DSDPCHKERR(info);*nzcols=0;}
        else if (row==m-1){
          info=DSDPVecZero(V);DSDPCHKERR(info);*nzcols=0;
          if (r){info=DSDPVecSetR(V,1.0);DSDPCHKERR(info);*nzcols=1;}
        } else if (M.dsdpops->matrownonzeros){
          info=DSDPVecGetSize(V,&m);DSDPCHKERR(info);
          info=DSDPVecGetArray(V,&cols);DSDPCHKERR(info);
          info=(M.dsdpops->matrownonzeros)(M.data,row-1,cols+1,nzcols,m-2); //DSDPChkMatError(M,info);
          info=DSDPVecRestoreArray(V,&cols);DSDPCHKERR(info);
          info=DSDPZeroFixedVariables(M,V);DSDPCHKERR(info);
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
            info=(A.dsdpops->matgetrank)(A.matdata,&rank,n); //DSDPChkDataError(A,info);
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
              info=(X.dsdpops->matzeroentries)(X.matdata); //DSDPChkMatError(X,info);
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
            info=SDPConeVecGetArray(V,&vv); DSDPCHKERR(info);
            info=SDPConeVecGetSize(V,&n); DSDPCHKERR(info);
            info=(A.dsdpops->matgeteig)(A.matdata,rr, eigenvalue, vv,n,S.indx+1,S.indx); //DSDPChkDataError(A,info);
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
            info=SDPConeVecGetSize(X,&n); DSDPCHKERR(info);
            info=SDPConeVecGetArray(B,&bb); DSDPCHKERR(info);
            info=SDPConeVecGetArray(X,&xx); DSDPCHKERR(info);
            info=(S.dsdpops->matinversemultiply)(S.matdata,IS.indx+1,IS.indx[0],bb,xx,n); //DSDPChkDMatError(S,info);
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
              info=SDPConeVecGetArray(V,&v); DSDPCHKERR(info);
              info=(X.dsdpops->mataddouterproduct)(X.matdata,alpha,v,n); //DSDPChkMatError(X,info);
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
                  info=SDPConeVecGetSize(W,&n); DSDPCHKERR(info);
                  info=SDPConeVecGetArray(W,&x); DSDPCHKERR(info);
                  info=(A.dsdpops->matvecvec)(A.matdata,x,n,v); // DSDPChkDataError(A,info);
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
              info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); // DSDPChkMatError(X,info);
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
              info=(X.dsdpops->matgetsize)(X.matdata,&n); //DSDPChkMatError(X,info);
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
              info=(X.dsdpops->matgeturarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
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
                  info=(blk[kk].ADATA.A[ii].dsdpops->matdot)(blk[kk].ADATA.A[ii].matdata,x,nn,n,&sum); 
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
              info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
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
	      info=(X.dsdpops->matscalediagonal)(X.matdata,dscale); // DSDPChkMatError(X,info);
	    } else {
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
          info=(M.dsdpops->mataddrow)(M.data,row-1,alpha,v+1,m-2); // DSDPChkMatError(M,info);
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
          if (S.dsdpops->matinversemultiply){
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
    	info=SDPConeVecDot(W,W2,&rtemp); DSDPCHKVARERR(ii,info);
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
            info=(X.dsdpops->matzeroentries)(X.matdata); //DSDPChkMatError(X,info);
          } else {
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
              info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
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
	      info=(X.dsdpops->matrestoreurarray)(X.matdata,v,&nn); //DSDPChkMatError(X,info);
	    } else {
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
	      } else {
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

