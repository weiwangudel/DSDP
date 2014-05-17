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
      info=DSDPSchurMatRowColumnScaling(M,i,Select,&ncols); DSDPCHKERR(info); 
      if (ncols==0){continue;}
   
      for (kt=0; kt<ATranspose.nnzblocks[i]; kt++){ /* Loop over all blocks */
        kk=ATranspose.nzblocks[i][kt];
        idA=ATranspose.idA[i][kt];
        info=DSDPBlockGetMatrix(&blk[kk].ADATA,idA,&ii,&scl,&AA);DSDPCHKBLOCKERR(kk,info);
        //Wei      if (ii!=i){DSDPSETERR2(8,"Data Transpose Error: var %d does not equal %d.\n",i,ii);}
        info = DSDPDataMatGetRank(AA,&rank,blk[kk].n);DSDPCHKBLOCKERR(kk,info);
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
        if (method1==DSDP_TRUE){info=DSDPVMatZeroEntries(T);DSDPCHKBLOCKERR(kk,info);}
        for (k=0; k<rank; k++){
  	
  	info=DSDPDataMatGetEig(AA,k,W,IS,&ack); DSDPCHKBLOCKERR(kk,info);
  	if (ack==0.0) continue;
  	ack*=scl;
  	info=DSDPDualMatInverseMultiply(S,IS,W,W2);DSDPCHKBLOCKERR(kk,info);
  
  	/* RHS terms */
  	info = SDPConeVecDot(W,W2,&rtemp); DSDPCHKBLOCKERR(kk,info);
  	if (rtemp==0.0) continue;
  	rhs1i+=rtemp*ack*bmu; rhs2i+=rtemp*ack*ggamma*mu;
  	ack*=(ggamma+bmu);
  
  	if (method1==DSDP_TRUE){
  	  info=DSDPVMatAddOuterProduct(T,ack*mu,W2);DSDPCHKBLOCKERR(kk,info);
  	} else {
  	  info=DSDPBlockvAv(&blk[kk].ADATA,ack*mu,Select,W2,MRowI);DSDPCHKBLOCKERR(kk,info);
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
          info=DSDPVMatScaleDiagonal(T,0.5); DSDPCHKERR(info);
          info=DSDPVMatGetSize(T, &n); DSDPCHKERR(info);
          info=DSDPVMatGetArray(T, &x, &nn); DSDPCHKERR(info);
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
          info=DSDPVMatRestoreArray(T, &x, &nn); DSDPCHKERR(info);
          info=DSDPVMatScaleDiagonal(T,2.0); DSDPCHKERR(info);
          //DSDPEventLogEnd(sdpdotevent);
          }   //Wei: end of inline function body of DSDPBlockADot
        }   /* End row computations for all of block ll     */
      }     /* End row computations for all blocks          */
      info=DSDPVecAddElement(vrhs1,i,rhs1i);DSDPCHKERR(info);
      info=DSDPVecAddElement(vrhs2,i,rhs2i);DSDPCHKERR(info);
      info=DSDPSchurMatAddRow(M,i,1.0,MRowI);DSDPCHKERR(info);
    }
  //}
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
    info=SDPConeComputeRHS(sdpcone,kk,mu,vrow,vrhs1,vrhs2); DSDPCHKBLOCKERR(kk,info);
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
      info=SDPConeComputeSS(sdpcone,kk,Y,blk[kk].T);DSDPCHKBLOCKERR(kk,info);
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

