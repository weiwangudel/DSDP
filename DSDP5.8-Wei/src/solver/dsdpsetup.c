#include "dsdp.h"
#include "dsdpsys.h"
#include "dsdp5.h"

#include "dsdpschurmat_impl.h"    //Inline need
#include "dsdpschurmat.h"
#include "dsdpbasictypes.h"
#include "dsdpsys.h"
#include "dsdpvec.h"
#include "dsdpbasictypes.h"
#include "dsdpcg.h"

//ComputeHessian Needed
#include "dsdpcone.h"    //Wei-inlining
#include "dsdpcone_impl.h"  //
#include <../sdp/dsdpsdp.h>    //Wei-inlining needed, rel



typedef enum { DSDPNoMatrix=1, DSDPUnfactoredMatrix=2, DSDPFactoredMatrix=3} DSDPCGType;
typedef struct{
  DSDPCGType type;
  DSDPSchurMat M;
  DSDPVec Diag;
  DSDP dsdp;
} DSDPCGMat;

//#include "libenergy.h"
/*!
  \file dsdpsetup.c
  \brief Create DSDP solver and its data strucutures.
*/

/*!
\fn int DSDPCreate(int m, DSDP *dsdpnew)
\brief Create a DSDP solver. FIRST DSDP routine!

\param m the number of variables y
\param *dsdpnew will be set to a new solver object
\sa DSDPSetDualObjective()
\sa DSDPSetup()
\sa DSDPSolve()
\sa DSDPDestroy()
\ingroup DSDPBasic

For example, to create a DSDP solver for a problem with 10 y variables,
\code
int m=10;
DSDP dsdp;
DSDPCreate(m,&dsdp);
\endcode
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPCreate"
int DSDPCreate(int m,DSDP* dsdpnew){

  DSDP dsdp;
  int info;

  DSDPFunctionBegin;

  DSDPCALLOC1(&dsdp,PD_DSDP,&info);DSDPCHKERR(info);
  *dsdpnew=dsdp;
  dsdp->keyid=DSDPKEY;

  /* Initialize some parameters */
  DSDPEventLogInitialize();
  dsdp->m=m;
  dsdp->maxcones=0;
  dsdp->ncones=0;
  dsdp->K=0;
  dsdp->setupcalled=DSDP_FALSE;
  dsdp->ybcone=0;
  dsdp->ndroutines=0;
  /*  info = DSDPSetStandardMonitor(dsdp);DSDPCHKERR(info);  */
  info = DSDPVecCreateSeq(m+2,&dsdp->b);DSDPCHKERR(info);
  info = DSDPVecZero(dsdp->b);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->b,&dsdp->y);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->b,&dsdp->ytemp);DSDPCHKERR(info);
  info = DSDPVecZero(dsdp->y);DSDPCHKERR(info);
  info = DSDPVecSetC(dsdp->y,-1.0);DSDPCHKERR(info);

  info = DSDPAddRCone(dsdp,&dsdp->rcone);DSDPCHKERR(info);
  info = DSDPCreateLUBoundsCone(dsdp,&dsdp->ybcone);DSDPCHKERR(info);

  info=DSDPSetDefaultStatistics(dsdp);DSDPCHKERR(info);
  info=DSDPSetDefaultParameters(dsdp);DSDPCHKERR(info);
  info=DSDPSetDefaultMonitors(dsdp);DSDPCHKERR(info);

  /*  info = DSDPMatInitialize(m,m,&dsdp->Q);DSDPCHKERR(info); */
  info = DSDPSchurMatInitialize(&dsdp->M);DSDPCHKERR(info);
  info = DSDPSetDefaultSchurMatrixStructure(dsdp); DSDPCHKERR(info);
  info = DSDPCGInitialize(&dsdp->sles); DSDPCHKERR(info);

  /* Set the one global variable
     sdat=dsdp;
  */
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDefaultStatistics"
/*! 
\fn int DSDPSetDefaultStatistics(DSDP dsdp);
\brief Set default statistics.
\param dsdp the solver
*/
int DSDPSetDefaultStatistics(DSDP dsdp){
  
  int i;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  dsdp->reason=CONTINUE_ITERATING;
  dsdp->pdfeasible=DSDP_PDUNKNOWN;
  dsdp->itnow=0;
  dsdp->pobj=  1.0e10;
  dsdp->ppobj=  1.0e10;
  dsdp->dobj= -1.0e+9;
  dsdp->ddobj= -1.0e+9;
  dsdp->dualitygap=dsdp->ppobj-dsdp->ddobj;
  dsdp->pstep=1.0;
  dsdp->dstep=0.0;
  for (i=0;i<MAX_XMAKERS;i++){
    dsdp->xmaker[i].mu=1.0e200;
    dsdp->xmaker[i].pstep=0.0;
  }
  dsdp->pnorm=0.001;
  dsdp->mu=1000.0;
  dsdp->np=0;
  dsdp->anorm=0;
  dsdp->bnorm=0;
  dsdp->cnorm=0;
  dsdp->tracex=0;
  dsdp->tracexs=0;
  dsdp->Mshift=0;
  dsdp->goty0=DSDP_FALSE;
  DSDPFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDefaultParameters"
/*! 
\fn int DSDPSetDefaultParameters(DSDP dsdp);
\brief Set default parameters.
\param dsdp the solver
*/
int DSDPSetDefaultParameters(DSDP dsdp){
  
  int info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);

  /* Stopping parameters */
  info=DSDPSetMaxIts(dsdp,500);DSDPCHKERR(info);
  info=DSDPSetGapTolerance(dsdp,1.0e-6);DSDPCHKERR(info);
  info=DSDPSetPNormTolerance(dsdp,1.0e30);DSDPCHKERR(info);
  if (dsdp->m<100){info=DSDPSetGapTolerance(dsdp,1.0e-7);DSDPCHKERR(info);}
  if (dsdp->m>3000){info=DSDPSetGapTolerance(dsdp,5.0e-6);DSDPCHKERR(info);}
  info=RConeSetType(dsdp->rcone,DSDPInfeasible);DSDPCHKERR(info);
  info=DSDPSetDualBound(dsdp,1.0e20);DSDPCHKERR(info);
  info=DSDPSetStepTolerance(dsdp,5.0e-2);DSDPCHKERR(info);
  info=DSDPSetRTolerance(dsdp,1.0e-6);DSDPCHKERR(info);
  info=DSDPSetPTolerance(dsdp,1.0e-4);DSDPCHKERR(info);
  /* Solver options */
  info=DSDPSetMaxTrustRadius(dsdp,1.0e10);DSDPCHKERR(info);
  info=DSDPUsePenalty(dsdp,0);DSDPCHKERR(info);
  info=DSDPSetInitialBarrierParameter(dsdp,-1.0);DSDPCHKERR(info);
  info=DSDPSetPotentialParameter(dsdp,3.0);DSDPCHKERR(info);
  info=DSDPUseDynamicRho(dsdp,1);DSDPCHKERR(info);
  info=DSDPSetR0(dsdp,-1.0);DSDPCHKERR(info);
  info=DSDPSetPenaltyParameter(dsdp,1.0e8);DSDPCHKERR(info);
  info=DSDPReuseMatrix(dsdp,4);DSDPCHKERR(info);
  if (dsdp->m>100){info=DSDPReuseMatrix(dsdp,7);DSDPCHKERR(info);}
  if (dsdp->m>1000){info=DSDPReuseMatrix(dsdp,10);DSDPCHKERR(info);}
  if (dsdp->m<=100){info=DSDPSetPotentialParameter(dsdp,5.0);DSDPCHKERR(info);DSDPCHKERR(info);}
  dsdp->maxschurshift=1.0e-11;
  dsdp->mu0=-1.0;
  dsdp->slestype=2;
  info = DSDPSetYBounds(dsdp,-1e7,1e7);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPSetDefaultMonitors"
/*! 
\fn int DSDPSetDefaultMonitors(DSDP dsdp);
\brief Set convergence monitor.
\param dsdp the solver
*/
int DSDPSetDefaultMonitors(DSDP dsdp){
  
  int info;

  DSDPFunctionBegin;
  DSDPValid(dsdp);
  dsdp->nmonitors=0;
  info=DSDPSetMonitor(dsdp,DSDPDefaultConvergence,(void*)&dsdp->conv); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetup(DSDP dsdp)
\brief Set up data structures in the solver and the cones associated with it.

\param dsdp the solver
\sa DSDPCreate()
\sa DSDPSolve()
\sa DSDPDestroy()
\ingroup DSDPBasic

This routine must be called before DSDPSolve().  Do not create
SDP, LP or other cones after calling this routines, and do not
set data into the cones after calling this routine.

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetUp"
int DSDPSetup(DSDP dsdp){
  
  int i,info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  
  /* Create the Work Vectors */
  info = DSDPVecDuplicate(dsdp->y,&dsdp->rhs1);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->rhs2);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->rhs);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->rhstemp);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->dy1);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->dy2);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->dy);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->y0);DSDPCHKERR(info);
  info = DSDPVecDuplicate(dsdp->y,&dsdp->xmakerrhs);DSDPCHKERR(info);
  for (i=0;i<MAX_XMAKERS;i++){
    info = DSDPVecDuplicate(dsdp->y,&dsdp->xmaker[i].y);DSDPCHKERR(info);
    info = DSDPVecDuplicate(dsdp->y,&dsdp->xmaker[i].dy);DSDPCHKERR(info);
    info = DSDPVecDuplicate(dsdp->y,&dsdp->xmaker[i].rhs);DSDPCHKERR(info);
  }

  /* Create M */
  info = DSDPSetUpCones(dsdp);DSDPCHKERR(info);
  info = DSDPSchurMatSetup(dsdp->M,dsdp->ytemp);DSDPCHKERR(info); 

  info = DSDPCGSetup(dsdp->sles,dsdp->ytemp); DSDPCHKERR(info);

  info = DSDPSetUpCones2(dsdp,dsdp->y,dsdp->M);DSDPCHKERR(info);
  info = DSDPGetConicDimension(dsdp,&dsdp->np);DSDPCHKERR(info);

  info=DSDPComputeDataNorms(dsdp);DSDPCHKERR(info);
  dsdp->pinfeas=dsdp->bnorm+1;
  dsdp->perror=dsdp->bnorm+1;
  info=DSDPScaleData(dsdp);DSDPCHKERR(info);

  info=DSDPGetConicDimension(dsdp,&dsdp->np);DSDPCHKERR(info);
  dsdp->solvetime=0;
  dsdp->cgtime=0;
  dsdp->ptime=0;
  dsdp->dtime=0;
  dsdp->ctime=0;
  info=DSDPEventLogRegister("Primal Step",&dsdp->ptime);
  info=DSDPEventLogRegister("Dual Step",&dsdp->dtime);
  info=DSDPEventLogRegister("Corrector Step",&dsdp->ctime);
  info=DSDPEventLogRegister("CG Solve",&dsdp->cgtime);
  info=DSDPEventLogRegister("DSDP Solve",&dsdp->solvetime);
  dsdp->setupcalled=DSDP_TRUE;
  DSDPFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DSDPGetSchurMatrix"
int DSDPGetSchurMatrix(DSDP dsdp, DSDPSchurMat *M){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *M=dsdp->M;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPGetConvergenceMonitor"
/*!
\fn int DSDPGetConvergenceMonitor(DSDP dsdp, ConvergenceMonitor**ctx);
\brief Get the structure containing convergence parameters.
\param dsdp the solver
\param *ctx will point to the structure.
\ingroup DSDPRoutines

\note 
This structure part of the DSDP structure.

 */
int DSDPGetConvergenceMonitor(DSDP dsdp, ConvergenceMonitor**ctx){
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  *ctx=&dsdp->conv;
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeDataNorms"
/*! 
\fn int DSDPComputeDataNorms(DSDP dsdp);
\brief Compute norms of A,C, and b.
\param dsdp the solver
*/
int DSDPComputeDataNorms(DSDP dsdp){
  int info;
  DSDPVec ytemp=dsdp->ytemp;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info = DSDPComputeANorm2(dsdp,ytemp);DSDPCHKERR(info);
  info = DSDPFixedVariablesNorm(dsdp->M,ytemp);DSDPCHKERR(info);
  info = DSDPVecGetC(ytemp,&dsdp->cnorm);DSDPCHKERR(info);
  dsdp->cnorm=sqrt(dsdp->cnorm);
  info = DSDPVecSetR(ytemp,0);DSDPCHKERR(info);
  info = DSDPVecSetC(ytemp,0);DSDPCHKERR(info);
  info = DSDPVecNorm1(ytemp,&dsdp->anorm);DSDPCHKERR(info);
  dsdp->anorm=sqrt(dsdp->anorm);
  DSDPLogInfo(0,2,"Norm of data: %4.2e\n",dsdp->anorm);
  info=DSDPVecCopy(dsdp->b,ytemp);DSDPCHKERR(info);
  info = DSDPVecSetR(ytemp,0);DSDPCHKERR(info);
  info = DSDPVecSetC(ytemp,0);DSDPCHKERR(info);
  info = DSDPVecNorm2(ytemp,&dsdp->bnorm);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPScaleData"
/*! 
\fn int DSDPScaleData(DSDP dsdp);
\brief Scale the matrix C.
\param dsdp the solver
*/
int DSDPScaleData(DSDP dsdp){
  int info;
  double scale;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  scale=1.0*dsdp->anorm;
  if (dsdp->bnorm){ scale/=dsdp->bnorm;}
  if (dsdp->cnorm){ scale/=dsdp->cnorm;}
  scale=DSDPMin(scale,1.0);
  scale=DSDPMax(scale,1.0e-6);
  if (dsdp->cnorm==0){  scale=1;}
  info=DSDPSetScale(dsdp,scale);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSolve(DSDP dsdp)
\brief Apply DSDP to the problem.

Call this routine after DSDPCreate() and DSDPSetup(),
and after setting the data.

\param dsdp is the solver
\sa DSDPCreate()
\sa DSDPGetSolutionType()
\sa DSDPGetDObjective()
\sa DSDPGetY()
\sa DSDPStopReason()
\ingroup DSDPBasic
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSolve"
int DSDPSolve(DSDP dsdp){
  int info;
  DSDPFunctionBegin;
  info=DSDPEventLogBegin(dsdp->solvetime);
  dsdp->pdfeasible=DSDP_PDUNKNOWN;
  info=DSDPSetConvergenceFlag(dsdp,CONTINUE_ITERATING);DSDPCHKERR(info); 
  info=DSDPInitializeVariables(dsdp);DSDPCHKERR(info);
  //info=DSDPSolveDynamicRho(dsdp);DSDPCHKERR(info);
  //Wei:   inlined function body of DSDPSolveDynamicRho 
  {
  //int DSDPSolveDynamicRho(DSDP dsdp){
  
    int info,attempt,maxattempts;
    double dd1,dd2,mutarget,ppnorm;
    DSDPTerminationReason reason;
    DSDPTruth cg1;
    DSDPTruth psdefinite;
  
    info=DSDPVecCopy(dsdp->y,dsdp->y0);DSDPCHKERR(info);
    for (dsdp->itnow=0; dsdp->itnow <= dsdp->maxiter+1 ; dsdp->itnow++){
  
      /* Check Convergence, and print information if desired */
      info=DSDPCheckConvergence(dsdp,&reason);DSDPCHKERR(info);
      if (reason != CONTINUE_ITERATING){break;}
      if (dsdp->mu0>0){dsdp->mutarget=DSDPMin(dsdp->mutarget,dsdp->mu0);}
      /* Compute the Gram matrix M and rhs */
      //info=DSDPComputeDualStepDirections(dsdp); DSDPCHKERR(info);
      //Wei: inlined function body of DSDPComputeDualStepDirections
      {
      //int DSDPComputeDualStepDirections(DSDP dsdp){
        int info,computem=1;
        double madd,ymax,cgtol=1e-7;
        DSDPTruth cg1,cg2,psdefinite;
      
        if (dsdp->itnow>30) dsdp->slestype=3;
        if (dsdp->rgap<1e-3) dsdp->slestype=3;
        if (dsdp->m<40) dsdp->slestype=3;
        if (0 && dsdp->itnow>20 && dsdp->m<500) dsdp->slestype=3;
        info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
        if (dsdp->slestype==1){
          cg1=DSDP_TRUE; cg2=DSDP_TRUE;
          info=DSDPInvertS(dsdp);DSDPCHKERR(info);
          info=DSDPComputeG(dsdp,dsdp->rhstemp,dsdp->rhs1,dsdp->rhs2);DSDPCHKERR(info);
          //info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs1,dsdp->dy1,cgtol,&cg1);DSDPCHKERR(info);
	  //Wei: inlined function body of DSDPCGSolve 
	  {
          //int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
	    DSDPSchurMat MM=dsdp->M;
	    DSDPVec RHS = dsdp->rhs1;
	    DSDPVec X = dsdp->dy1;
            DSDPTruth *success = &cg1;

            int iter=0,n,info,maxit=10;
            double dd,ymax;
            DSDPCG *sles=dsdp->sles; 
            DSDPCGMat    CGM;
          
            info=DSDPEventLogBegin(dsdp->cgtime);
            info=DSDPVecZero(X);DSDPCHKERR(info);
            info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
            *success=DSDP_TRUE;
            if (0){
              maxit=0;
            } else if (dsdp->slestype==1){
          
              CGM.type=DSDPNoMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              cgtol*=1000;
              maxit=5;
          
            } else if (dsdp->slestype==2){
              CGM.type=DSDPUnfactoredMatrix;
              CGM.M=MM;
              CGM.Diag=sles->Diag;
              CGM.dsdp=dsdp;
              cgtol*=100;
              maxit=(int)sqrt(1.0*n)+10;
              if (maxit>20) maxit=20;
              info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
              /*
                info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
                info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
              */
              
            } else if (dsdp->slestype==3){
              CGM.type=DSDPFactoredMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              maxit=0;
              info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
              if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
              if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
                maxit=6;
              } else if (dsdp->rgap<1e-5){
                maxit=3;
              }
          
              //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
              //Wei: inlined function body of DSDPSchurMatSolve
              {
                
              //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                int info;
                //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                //Wei: inline function body DSDPSchurMatSolveM
                {  
                //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                //{
                  int info,n;
                  double *xx,*bb;
                  //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                  if (MM.dsdpops->matsolve){
                    info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                    info=DSDPVecZero(X);DSDPCHKERR(info);
                    info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                    info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                    info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                  } else {
                    //DSDPNoOperationError(MM);
          	  exit(-1);         // omit error handling
                  }
                  info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                  info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                  //info=DSDPEventLogEnd(hsolveevent);
                //}
                }   // end of inlined function body DSDPSchurMatSolveM
                info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
              //} // end of original function body
          
              } // end of inlined function DSDPSchurMatSolve
          
            } else if (dsdp->slestype==4){
              CGM.type=DSDPFactoredMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              maxit=3;
              //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
              //Wei: inlined function body of DSDPSchurMatSolve
              {
                
              //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                int info;
                //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                //Wei: inline function body DSDPSchurMatSolveM
                {  
                //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                //{
                  int info,n;
                  double *xx,*bb;
                  //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                  if (MM.dsdpops->matsolve){
                    info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                    info=DSDPVecZero(X);DSDPCHKERR(info);
                    info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                    info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                    info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                  } else {
                    //DSDPNoOperationError(MM);
          	  exit(-1);         // omit error handling
                  }
                  info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                  info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                  //info=DSDPEventLogEnd(hsolveevent);
                //}
                }   // end of inlined function body DSDPSchurMatSolveM
                info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
              //} // end of original function body
          
              } // end of inlined function DSDPSchurMatSolve
            }
            if (n<maxit) maxit=n;
            
            info=DSDPConjugateGradient(CGM,X,RHS,
          			     sles->R,sles->BR,sles->P,sles->BP,
          			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);
          
            if (iter>=maxit) *success=DSDP_FALSE;
            info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
            if (dd<0) *success=DSDP_FALSE;
            info=DSDPEventLogEnd(dsdp->cgtime);
          //}


          } // end of inlined function body of DSDPCGSolve
          if (cg1==DSDP_TRUE){
	    //info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs2,dsdp->dy2,cgtol,&cg2);DSDPCHKERR(info);
	    {
            //int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
	      DSDPSchurMat MM=dsdp->M;
	      DSDPVec RHS = dsdp->rhs2;
	      DSDPVec X = dsdp->dy2;
              DSDPTruth *success = &cg2;

              int iter=0,n,info,maxit=10;
              double dd,ymax;
              DSDPCG *sles=dsdp->sles; 
              DSDPCGMat    CGM;
            
              info=DSDPEventLogBegin(dsdp->cgtime);
              info=DSDPVecZero(X);DSDPCHKERR(info);
              info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
              *success=DSDP_TRUE;
              if (0){
                maxit=0;
              } else if (dsdp->slestype==1){
            
                CGM.type=DSDPNoMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                cgtol*=1000;
                maxit=5;
            
              } else if (dsdp->slestype==2){
                CGM.type=DSDPUnfactoredMatrix;
                CGM.M=MM;
                CGM.Diag=sles->Diag;
                CGM.dsdp=dsdp;
                cgtol*=100;
                maxit=(int)sqrt(1.0*n)+10;
                if (maxit>20) maxit=20;
                info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
                /*
                  info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
                  info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
                */
                
              } else if (dsdp->slestype==3){
                CGM.type=DSDPFactoredMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                maxit=0;
                info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
                if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
                if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
                  maxit=6;
                } else if (dsdp->rgap<1e-5){
                  maxit=3;
                }
            
                //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
                //Wei: inlined function body of DSDPSchurMatSolve
                {
                  
                //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                  int info;
                  //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                  //Wei: inline function body DSDPSchurMatSolveM
                  {  
                  //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                  //{
                    int info,n;
                    double *xx,*bb;
                    //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                    if (MM.dsdpops->matsolve){
                      info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                      info=DSDPVecZero(X);DSDPCHKERR(info);
                      info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                      info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                      info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                    } else {
                      //DSDPNoOperationError(MM);
            	  exit(-1);         // omit error handling
                    }
                    info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                    info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                    //info=DSDPEventLogEnd(hsolveevent);
                  //}
                  }   // end of inlined function body DSDPSchurMatSolveM
                  info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                  info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
                //} // end of original function body
            
                } // end of inlined function DSDPSchurMatSolve
            
              } else if (dsdp->slestype==4){
                CGM.type=DSDPFactoredMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                maxit=3;
                //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
                //Wei: inlined function body of DSDPSchurMatSolve
                {
                  
                //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                  int info;
                  //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                  //Wei: inline function body DSDPSchurMatSolveM
                  {  
                  //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                  //{
                    int info,n;
                    double *xx,*bb;
                    //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                    if (MM.dsdpops->matsolve){
                      info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                      info=DSDPVecZero(X);DSDPCHKERR(info);
                      info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                      info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                      info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                    } else {
                      //DSDPNoOperationError(MM);
            	  exit(-1);         // omit error handling
                    }
                    info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                    info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                    //info=DSDPEventLogEnd(hsolveevent);
                  //}
                  }   // end of inlined function body DSDPSchurMatSolveM
                  info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                  info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
                //} // end of original function body
            
                } // end of inlined function DSDPSchurMatSolve
              }
              if (n<maxit) maxit=n;
              
              info=DSDPConjugateGradient(CGM,X,RHS,
            			     sles->R,sles->BR,sles->P,sles->BP,
            			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);
            
              if (iter>=maxit) *success=DSDP_FALSE;
              info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
              if (dd<0) *success=DSDP_FALSE;
              info=DSDPEventLogEnd(dsdp->cgtime);
            //}


            } // end of inlined function body of DSDPCGSolve
	  }
          if (cg1==DSDP_FALSE || cg2==DSDP_FALSE) dsdp->slestype=2;
        }
        if (dsdp->slestype==2){
          cg1=DSDP_TRUE; cg2=DSDP_TRUE;
          DSDPLogInfo(0,9,"Compute Hessian\n");
          info=DSDPInvertS(dsdp);DSDPCHKERR(info);
          //info=DSDPComputeHessian(dsdp,dsdp->M,dsdp->rhs1,dsdp->rhs2);DSDPCHKERR(info);
          //Wei: inline function body of DSDPComputeHessian
          {
          //int DSDPComputeHessian( DSDP dsdp , DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2)
          //{
            int info,kk; double r;
            //DSDPEventLogBegin(ConeComputeH); // ignore log, do not know how to handle staic ConeComputeH
            dsdp->schurmu=dsdp->mutarget;
            info=DSDPVecGetR(dsdp->y,&r);DSDPCHKERR(info);
            info=DSDPSchurMatSetR(dsdp->M,r);DSDPCHKERR(info);
            info=DSDPSchurMatZeroEntries(dsdp->M);DSDPCHKERR(info);
            info=DSDPVecZero(dsdp->rhs1);DSDPCHKERR(info);
            info=DSDPVecZero(dsdp->rhs2);DSDPCHKERR(info);
            info=DSDPVecZero((dsdp->M).schur->rhs3);DSDPCHKERR(info);
            info=DSDPObjectiveGH(dsdp,dsdp->M,dsdp->rhs1); DSDPCHKERR(info);
            for (kk=dsdp->ncones-1;kk>=0;kk--){
              DSDPEventLogBegin(dsdp->K[kk].coneid);
              //info=DSDPConeComputeHessian(dsdp->K[kk].cone,dsdp->schurmu,M,vrhs1,vrhs2);DSDPCHKCONEERR(kk,info);
              //Wei: Inline(1) Function Body of DSDPConeComputeHessian
              if (dsdp->K[kk].cone.dsdpops->conehessian){
                // Wei@05/15/14: printf("%x\n", dsdp->K[kk].cone.dsdpops->conehessian);
                // Wei: Cannot inline: conehessian are dynamically changed. 
                // Wei: KSDPConeComputeHessian only takes 1/3 of calls, but consumes most of time 
                dsdp->K[kk].cone.dsdpops->conehessian(dsdp->K[kk].cone.conedata,dsdp->schurmu,dsdp->M,dsdp->rhs1,dsdp->rhs2); // DSDPChkConeError(K,info);
              } else {
                 //DSDPNoOperationError(K);
                 // error handling omitted
                 exit(-1);
              }
              DSDPEventLogEnd(dsdp->K[kk].coneid);
            }
            info=DSDPSchurMatAssemble(dsdp->M);DSDPCHKERR(info);
            /*    DSDPSchurMatView(M); */
            info=DSDPSchurMatReducePVec(dsdp->M,dsdp->rhs1);DSDPCHKERR(info);
            info=DSDPSchurMatReducePVec(dsdp->M,dsdp->rhs2);DSDPCHKERR(info);
            info=DSDPSchurMatReducePVec(dsdp->M,(dsdp->M).schur->rhs3);DSDPCHKERR(info);
            if (0 && dsdp->UsePenalty==DSDPNever){   //why having this 0 && ????
              info=DSDPVecAXPY(1.0,dsdp->M.schur->rhs3,dsdp->rhs2);DSDPCHKERR(info);
              info=DSDPVecZero((dsdp->M).schur->rhs3);DSDPCHKERR(info);
              info=DSDPVecZero((dsdp->M).schur->dy3);DSDPCHKERR(info);
              info=DSDPVecSetR(dsdp->rhs1,0);DSDPCHKERR(info);
              info=DSDPVecSetR(dsdp->rhs2,r);DSDPCHKERR(info);
            }
            //DSDPEventLogEnd(ConeComputeH);
          //}
      
          } // end of inline function body of DSDPComputeHessian
          computem=0;
          DSDPLogInfo(0,9,"Apply CG\n");
          //info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs1,dsdp->dy1,cgtol,&cg1);DSDPCHKERR(info);
	  {
          //int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
	    DSDPSchurMat MM=dsdp->M;
	    DSDPVec RHS = dsdp->rhs1;
	    DSDPVec X = dsdp->dy1;
            DSDPTruth *success = &cg1;

            int iter=0,n,info,maxit=10;
            double dd,ymax;
            DSDPCG *sles=dsdp->sles; 
            DSDPCGMat    CGM;
          
            info=DSDPEventLogBegin(dsdp->cgtime);
            info=DSDPVecZero(X);DSDPCHKERR(info);
            info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
            *success=DSDP_TRUE;
            if (0){
              maxit=0;
            } else if (dsdp->slestype==1){
          
              CGM.type=DSDPNoMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              cgtol*=1000;
              maxit=5;
          
            } else if (dsdp->slestype==2){
              CGM.type=DSDPUnfactoredMatrix;
              CGM.M=MM;
              CGM.Diag=sles->Diag;
              CGM.dsdp=dsdp;
              cgtol*=100;
              maxit=(int)sqrt(1.0*n)+10;
              if (maxit>20) maxit=20;
              info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
              /*
                info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
                info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
              */
              
            } else if (dsdp->slestype==3){
              CGM.type=DSDPFactoredMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              maxit=0;
              info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
              if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
              if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
                maxit=6;
              } else if (dsdp->rgap<1e-5){
                maxit=3;
              }
          
              //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
              //Wei: inlined function body of DSDPSchurMatSolve
              {
                
              //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                int info;
                //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                //Wei: inline function body DSDPSchurMatSolveM
                {  
                //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                //{
                  int info,n;
                  double *xx,*bb;
                  //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                  if (MM.dsdpops->matsolve){
                    info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                    info=DSDPVecZero(X);DSDPCHKERR(info);
                    info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                    info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                    info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                  } else {
                    //DSDPNoOperationError(MM);
          	  exit(-1);         // omit error handling
                  }
                  info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                  info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                  //info=DSDPEventLogEnd(hsolveevent);
                //}
                }   // end of inlined function body DSDPSchurMatSolveM
                info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
              //} // end of original function body
          
              } // end of inlined function DSDPSchurMatSolve
          
            } else if (dsdp->slestype==4){
              CGM.type=DSDPFactoredMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              maxit=3;
              //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
              //Wei: inlined function body of DSDPSchurMatSolve
              {
                
              //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                int info;
                //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                //Wei: inline function body DSDPSchurMatSolveM
                {  
                //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                //{
                  int info,n;
                  double *xx,*bb;
                  //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                  if (MM.dsdpops->matsolve){
                    info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                    info=DSDPVecZero(X);DSDPCHKERR(info);
                    info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                    info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                    info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                  } else {
                    //DSDPNoOperationError(MM);
          	  exit(-1);         // omit error handling
                  }
                  info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                  info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                  //info=DSDPEventLogEnd(hsolveevent);
                //}
                }   // end of inlined function body DSDPSchurMatSolveM
                info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
              //} // end of original function body
          
              } // end of inlined function DSDPSchurMatSolve
            }
            if (n<maxit) maxit=n;
            
            info=DSDPConjugateGradient(CGM,X,RHS,
          			     sles->R,sles->BR,sles->P,sles->BP,
          			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);
          
            if (iter>=maxit) *success=DSDP_FALSE;
            info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
            if (dd<0) *success=DSDP_FALSE;
            info=DSDPEventLogEnd(dsdp->cgtime);
          //}


          } // end of inlined function body of DSDPCGSolve
          if (cg1==DSDP_TRUE) {
	    //info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs2,dsdp->dy2,cgtol,&cg2);DSDPCHKERR(info);
	    {
            //int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
	      DSDPSchurMat MM=dsdp->M;
	      DSDPVec RHS = dsdp->rhs2;
	      DSDPVec X = dsdp->dy2;
              DSDPTruth *success = &cg2;

              int iter=0,n,info,maxit=10;
              double dd,ymax;
              DSDPCG *sles=dsdp->sles; 
              DSDPCGMat    CGM;
            
              info=DSDPEventLogBegin(dsdp->cgtime);
              info=DSDPVecZero(X);DSDPCHKERR(info);
              info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
              *success=DSDP_TRUE;
              if (0){
                maxit=0;
              } else if (dsdp->slestype==1){
            
                CGM.type=DSDPNoMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                cgtol*=1000;
                maxit=5;
            
              } else if (dsdp->slestype==2){
                CGM.type=DSDPUnfactoredMatrix;
                CGM.M=MM;
                CGM.Diag=sles->Diag;
                CGM.dsdp=dsdp;
                cgtol*=100;
                maxit=(int)sqrt(1.0*n)+10;
                if (maxit>20) maxit=20;
                info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
                /*
                  info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
                  info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
                */
                
              } else if (dsdp->slestype==3){
                CGM.type=DSDPFactoredMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                maxit=0;
                info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
                if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
                if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
                  maxit=6;
                } else if (dsdp->rgap<1e-5){
                  maxit=3;
                }
            
                //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
                //Wei: inlined function body of DSDPSchurMatSolve
                {
                  
                //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                  int info;
                  //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                  //Wei: inline function body DSDPSchurMatSolveM
                  {  
                  //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                  //{
                    int info,n;
                    double *xx,*bb;
                    //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                    if (MM.dsdpops->matsolve){
                      info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                      info=DSDPVecZero(X);DSDPCHKERR(info);
                      info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                      info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                      info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                    } else {
                      //DSDPNoOperationError(MM);
            	  exit(-1);         // omit error handling
                    }
                    info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                    info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                    //info=DSDPEventLogEnd(hsolveevent);
                  //}
                  }   // end of inlined function body DSDPSchurMatSolveM
                  info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                  info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
                //} // end of original function body
            
                } // end of inlined function DSDPSchurMatSolve
            
              } else if (dsdp->slestype==4){
                CGM.type=DSDPFactoredMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                maxit=3;
                //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
                //Wei: inlined function body of DSDPSchurMatSolve
                {
                  
                //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                  int info;
                  //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                  //Wei: inline function body DSDPSchurMatSolveM
                  {  
                  //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                  //{
                    int info,n;
                    double *xx,*bb;
                    //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                    if (MM.dsdpops->matsolve){
                      info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                      info=DSDPVecZero(X);DSDPCHKERR(info);
                      info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                      info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                      info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                    } else {
                      //DSDPNoOperationError(MM);
            	  exit(-1);         // omit error handling
                    }
                    info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                    info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                    //info=DSDPEventLogEnd(hsolveevent);
                  //}
                  }   // end of inlined function body DSDPSchurMatSolveM
                  info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                  info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
                //} // end of original function body
            
                } // end of inlined function DSDPSchurMatSolve
              }
              if (n<maxit) maxit=n;
              
              info=DSDPConjugateGradient(CGM,X,RHS,
            			     sles->R,sles->BR,sles->P,sles->BP,
            			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);
            
              if (iter>=maxit) *success=DSDP_FALSE;
              info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
              if (dd<0) *success=DSDP_FALSE;
              info=DSDPEventLogEnd(dsdp->cgtime);
            //}


            } // end of inlined function body of DSDPCGSolve
	  }
          if (cg1==DSDP_FALSE || cg2==DSDP_FALSE) dsdp->slestype=3;
          
        }
        if (dsdp->slestype==3){
          DSDPLogInfo(0,9,"Factor Hessian\n");
          psdefinite=DSDP_FALSE;
          if (dsdp->Mshift < 1e-12 || dsdp->rgap<0.1 || dsdp->Mshift > 1e-6){
            madd=dsdp->Mshift;
          } else {
            madd=1e-13;
          }
          if (computem){
            info=DSDPInvertS(dsdp);DSDPCHKERR(info);
          }
          while (psdefinite==DSDP_FALSE){
            if (0==1 && dsdp->Mshift>dsdp->maxschurshift){ 
      	info = DSDPSetConvergenceFlag(dsdp,DSDP_INDEFINITE_SCHUR_MATRIX); DSDPCHKERR(info);  
      	break;
            }
            if (0 && dsdp->Mshift*ymax>dsdp->pinfeastol/10){ 
      	info = DSDPSetConvergenceFlag(dsdp,DSDP_INDEFINITE_SCHUR_MATRIX); DSDPCHKERR(info);  
      	break;
            }
            if (madd*ymax>dsdp->pinfeastol*1000){ 
      	info = DSDPSetConvergenceFlag(dsdp,DSDP_INDEFINITE_SCHUR_MATRIX); DSDPCHKERR(info);  
      	break;
            }
            if (computem){
      	//info=DSDPComputeHessian(dsdp,dsdp->M,dsdp->rhs1,dsdp->rhs2);DSDPCHKERR(info);
              //Wei: inline function body of DSDPComputeHessian
              {
              //int DSDPComputeHessian( DSDP dsdp , DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2)
              //{
                int info,kk; double r;
                //DSDPEventLogBegin(ConeComputeH); // ignore log, do not know how to handle staic ConeComputeH
                dsdp->schurmu=dsdp->mutarget;
                info=DSDPVecGetR(dsdp->y,&r);DSDPCHKERR(info);
                info=DSDPSchurMatSetR(dsdp->M,r);DSDPCHKERR(info);
                info=DSDPSchurMatZeroEntries(dsdp->M);DSDPCHKERR(info);
                info=DSDPVecZero(dsdp->rhs1);DSDPCHKERR(info);
                info=DSDPVecZero(dsdp->rhs2);DSDPCHKERR(info);
                info=DSDPVecZero((dsdp->M).schur->rhs3);DSDPCHKERR(info);
                info=DSDPObjectiveGH(dsdp,dsdp->M,dsdp->rhs1); DSDPCHKERR(info);
                for (kk=dsdp->ncones-1;kk>=0;kk--){
                  DSDPEventLogBegin(dsdp->K[kk].coneid);
                  //info=DSDPConeComputeHessian(dsdp->K[kk].cone,dsdp->schurmu,M,vrhs1,vrhs2);DSDPCHKCONEERR(kk,info);
                  //Wei: Inline(1) Function Body of DSDPConeComputeHessian
                  if (dsdp->K[kk].cone.dsdpops->conehessian){
                    // Wei@05/15/14: printf("%x\n", dsdp->K[kk].cone.dsdpops->conehessian);
                    // Wei: Cannot inline: conehessian are dynamically changed. 
                    // Wei: KSDPConeComputeHessian only takes 1/3 of calls, but consumes most of time 
                    dsdp->K[kk].cone.dsdpops->conehessian(dsdp->K[kk].cone.conedata,dsdp->schurmu,dsdp->M,dsdp->rhs1,dsdp->rhs2); // DSDPChkConeError(K,info);
                  } else {
                     //DSDPNoOperationError(K);
                     // error handling omitted
                     exit(-1);
                  }
                  DSDPEventLogEnd(dsdp->K[kk].coneid);
                }
                info=DSDPSchurMatAssemble(dsdp->M);DSDPCHKERR(info);
                /*    DSDPSchurMatView(M); */
                info=DSDPSchurMatReducePVec(dsdp->M,dsdp->rhs1);DSDPCHKERR(info);
                info=DSDPSchurMatReducePVec(dsdp->M,dsdp->rhs2);DSDPCHKERR(info);
                info=DSDPSchurMatReducePVec(dsdp->M,(dsdp->M).schur->rhs3);DSDPCHKERR(info);
                if (0 && dsdp->UsePenalty==DSDPNever){   //why having this 0 && ????
                  info=DSDPVecAXPY(1.0,dsdp->M.schur->rhs3,dsdp->rhs2);DSDPCHKERR(info);
                  info=DSDPVecZero((dsdp->M).schur->rhs3);DSDPCHKERR(info);
                  info=DSDPVecZero((dsdp->M).schur->dy3);DSDPCHKERR(info);
                  info=DSDPVecSetR(dsdp->rhs1,0);DSDPCHKERR(info);
                  info=DSDPVecSetR(dsdp->rhs2,r);DSDPCHKERR(info);
                }
                //DSDPEventLogEnd(ConeComputeH);
              //}
      
          	} // end of inline function body of DSDPComputeHessian
      	
            }
            if (0==1){info=DSDPSchurMatView(dsdp->M);DSDPCHKERR(info);}
            info = DSDPSchurMatShiftDiagonal(dsdp->M,madd);DSDPCHKERR(info);
            //info = DSDPSchurMatFactor(dsdp->M,&psdefinite); DSDPCHKERR(info);
            //Wei: inline function body of DSDPSchurMatFactor 
            {
              //int DSDPSchurMatFactor(DSDPSchurMat M, DSDPTruth *successful)
              //{
              int info,flag=0;
              DSDPVec rhs3=(dsdp->M).schur->rhs3,dy3=(dsdp->M).schur->dy3;
      
              psdefinite=DSDP_TRUE;
              //DSDPEventLogBegin(hfactorevent);  //cannot link hfactorevent
              if ((dsdp->M).dsdpops->matfactor){
                info=((dsdp->M).dsdpops->matfactor)((dsdp->M).data,&flag); // DSDPChkMatError((dsdp->M),info);
                if (flag){ 
                  psdefinite=DSDP_FALSE;
                  DSDPLogInfo(0,2,"Indefinite Schur Matrix -- Bad Factorization\n");
                }
              } else {
                //DSDPNoOperationError(M);
                exit(-1);   // omit error handling
              }
              //DSDPEventLogEnd(hfactorevent);
              if ((dsdp->M).schur->r){
                info=DSDPSchurMatSolveM((dsdp->M),rhs3,dy3);DSDPCHKERR(info);
              }
              else {info=DSDPVecZero(dy3);DSDPCHKERR(info);}
              //}
            }  // end of inline function body of DSDPSchurMatFactor 
            computem=1;
            if (psdefinite==DSDP_FALSE){ 
      	madd=madd*4 + 1.0e-13;
            }
          }
          dsdp->Mshift=madd;
          if (psdefinite==DSDP_TRUE){
          //info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs1,dsdp->dy1,cgtol,&cg1);DSDPCHKERR(info);
	  {
          //int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
	    DSDPSchurMat MM=dsdp->M;
	    DSDPVec RHS = dsdp->rhs1;
	    DSDPVec X = dsdp->dy1;
            DSDPTruth *success = &cg1;

            int iter=0,n,info,maxit=10;
            double dd,ymax;
            DSDPCG *sles=dsdp->sles; 
            DSDPCGMat    CGM;
          
            info=DSDPEventLogBegin(dsdp->cgtime);
            info=DSDPVecZero(X);DSDPCHKERR(info);
            info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
            *success=DSDP_TRUE;
            if (0){
              maxit=0;
            } else if (dsdp->slestype==1){
          
              CGM.type=DSDPNoMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              cgtol*=1000;
              maxit=5;
          
            } else if (dsdp->slestype==2){
              CGM.type=DSDPUnfactoredMatrix;
              CGM.M=MM;
              CGM.Diag=sles->Diag;
              CGM.dsdp=dsdp;
              cgtol*=100;
              maxit=(int)sqrt(1.0*n)+10;
              if (maxit>20) maxit=20;
              info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
              /*
                info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
                info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
              */
              
            } else if (dsdp->slestype==3){
              CGM.type=DSDPFactoredMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              maxit=0;
              info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
              if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
              if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
                maxit=6;
              } else if (dsdp->rgap<1e-5){
                maxit=3;
              }
          
              //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
              //Wei: inlined function body of DSDPSchurMatSolve
              {
                
              //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                int info;
                //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                //Wei: inline function body DSDPSchurMatSolveM
                {  
                //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                //{
                  int info,n;
                  double *xx,*bb;
                  //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                  if (MM.dsdpops->matsolve){
                    info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                    info=DSDPVecZero(X);DSDPCHKERR(info);
                    info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                    info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                    info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                  } else {
                    //DSDPNoOperationError(MM);
          	  exit(-1);         // omit error handling
                  }
                  info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                  info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                  //info=DSDPEventLogEnd(hsolveevent);
                //}
                }   // end of inlined function body DSDPSchurMatSolveM
                info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
              //} // end of original function body
          
              } // end of inlined function DSDPSchurMatSolve
          
            } else if (dsdp->slestype==4){
              CGM.type=DSDPFactoredMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              maxit=3;
              //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
              //Wei: inlined function body of DSDPSchurMatSolve
              {
                
              //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                int info;
                //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                //Wei: inline function body DSDPSchurMatSolveM
                {  
                //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                //{
                  int info,n;
                  double *xx,*bb;
                  //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                  if (MM.dsdpops->matsolve){
                    info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                    info=DSDPVecZero(X);DSDPCHKERR(info);
                    info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                    info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                    info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                  } else {
                    //DSDPNoOperationError(MM);
          	  exit(-1);         // omit error handling
                  }
                  info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                  info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                  //info=DSDPEventLogEnd(hsolveevent);
                //}
                }   // end of inlined function body DSDPSchurMatSolveM
                info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
              //} // end of original function body
          
              } // end of inlined function DSDPSchurMatSolve
            }
            if (n<maxit) maxit=n;
            
            info=DSDPConjugateGradient(CGM,X,RHS,
          			     sles->R,sles->BR,sles->P,sles->BP,
          			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);
          
            if (iter>=maxit) *success=DSDP_FALSE;
            info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
            if (dd<0) *success=DSDP_FALSE;
            info=DSDPEventLogEnd(dsdp->cgtime);
          //}


          } // end of inlined function body of DSDPCGSolve
	    
            //info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs2,dsdp->dy2,cgtol,&cg2);DSDPCHKERR(info);
	    {
            //int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
	      DSDPSchurMat MM=dsdp->M;
	      DSDPVec RHS = dsdp->rhs2;
	      DSDPVec X = dsdp->dy2;
              DSDPTruth *success = &cg2;

              int iter=0,n,info,maxit=10;
              double dd,ymax;
              DSDPCG *sles=dsdp->sles; 
              DSDPCGMat    CGM;
            
              info=DSDPEventLogBegin(dsdp->cgtime);
              info=DSDPVecZero(X);DSDPCHKERR(info);
              info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
              *success=DSDP_TRUE;
              if (0){
                maxit=0;
              } else if (dsdp->slestype==1){
            
                CGM.type=DSDPNoMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                cgtol*=1000;
                maxit=5;
            
              } else if (dsdp->slestype==2){
                CGM.type=DSDPUnfactoredMatrix;
                CGM.M=MM;
                CGM.Diag=sles->Diag;
                CGM.dsdp=dsdp;
                cgtol*=100;
                maxit=(int)sqrt(1.0*n)+10;
                if (maxit>20) maxit=20;
                info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
                /*
                  info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
                  info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
                */
                
              } else if (dsdp->slestype==3){
                CGM.type=DSDPFactoredMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                maxit=0;
                info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
                if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
                if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
                  maxit=6;
                } else if (dsdp->rgap<1e-5){
                  maxit=3;
                }
            
                //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
                //Wei: inlined function body of DSDPSchurMatSolve
                {
                  
                //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                  int info;
                  //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                  //Wei: inline function body DSDPSchurMatSolveM
                  {  
                  //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                  //{
                    int info,n;
                    double *xx,*bb;
                    //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                    if (MM.dsdpops->matsolve){
                      info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                      info=DSDPVecZero(X);DSDPCHKERR(info);
                      info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                      info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                      info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                    } else {
                      //DSDPNoOperationError(MM);
            	  exit(-1);         // omit error handling
                    }
                    info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                    info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                    //info=DSDPEventLogEnd(hsolveevent);
                  //}
                  }   // end of inlined function body DSDPSchurMatSolveM
                  info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                  info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
                //} // end of original function body
            
                } // end of inlined function DSDPSchurMatSolve
            
              } else if (dsdp->slestype==4){
                CGM.type=DSDPFactoredMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                maxit=3;
                //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
                //Wei: inlined function body of DSDPSchurMatSolve
                {
                  
                //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                  int info;
                  //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                  //Wei: inline function body DSDPSchurMatSolveM
                  {  
                  //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                  //{
                    int info,n;
                    double *xx,*bb;
                    //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                    if (MM.dsdpops->matsolve){
                      info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                      info=DSDPVecZero(X);DSDPCHKERR(info);
                      info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                      info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                      info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                    } else {
                      //DSDPNoOperationError(MM);
            	  exit(-1);         // omit error handling
                    }
                    info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                    info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                    //info=DSDPEventLogEnd(hsolveevent);
                  //}
                  }   // end of inlined function body DSDPSchurMatSolveM
                  info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                  info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
                //} // end of original function body
            
                } // end of inlined function DSDPSchurMatSolve
              }
              if (n<maxit) maxit=n;
              
              info=DSDPConjugateGradient(CGM,X,RHS,
            			     sles->R,sles->BR,sles->P,sles->BP,
            			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);
            
              if (iter>=maxit) *success=DSDP_FALSE;
              info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
              if (dd<0) *success=DSDP_FALSE;
              info=DSDPEventLogEnd(dsdp->cgtime);
            //}


            } // end of inlined function body of DSDPCGSolve
	    
          }
        }
        
      //}

      }  // end of inlined function body of DSDPComputeDualStepDirections
      if (dsdp->reason==DSDP_INDEFINITE_SCHUR_MATRIX){continue;}
  
      info=DSDPComputePDY(dsdp,dsdp->mutarget,dsdp->dy,&dsdp->pnorm); DSDPCHKERR(info);
  
      DSDPEventLogBegin(dsdp->ptime);
      info=DSDPComputePY(dsdp,1.0,dsdp->ytemp);DSDPCHKERR(info);
      //info=DSDPComputeSS(dsdp,dsdp->ytemp,PRIMAL_FACTOR,&psdefinite);DSDPCHKERR(info);
      //Wei: inlined function body of DSDPComputeSS 
      {
      //int DSDPComputeSS(DSDP dsdp, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
        DSDPVec Y=dsdp->ytemp;
        DSDPDualFactorMatrix flag = PRIMAL_FACTOR;
        DSDPTruth *ispsdefinite = &psdefinite;

        int info,kk;
        DSDPTruth psd=DSDP_TRUE;
        if (flag==DUAL_FACTOR){
          //DSDPEventLogBegin(ConeComputeS);
        } else if (flag==PRIMAL_FACTOR){
          //DSDPEventLogBegin(ConeComputeSS);
        }
        for (kk=dsdp->ncones-1; kk>=0 && psd==DSDP_TRUE;kk--){
          DSDPEventLogBegin(dsdp->K[kk].coneid);
          //info=DSDPConeComputeS(dsdp->K[kk].cone,Y,flag,&psd); DSDPCHKCONEERR(kk,info);
          //Wei: inlined function body of DSDPConeComputeS
          {
          //int DSDPConeComputeS(DSDPCone K, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
            DSDPCone K=dsdp->K[kk].cone;
            DSDPTruth *ispsdefinite = &psd;
            int info;
            if (K.dsdpops->conecomputes){
              info=K.dsdpops->conecomputes(K.conedata,Y,flag,ispsdefinite); //DSDPChkConeError(K,info);
            } else {
              //DSDPNoOperationError(K);
      	exit(-1);   //omit error handling
            }
          //}
          } // end of inlined function body of DSDPConeComputeS
          DSDPEventLogEnd(dsdp->K[kk].coneid);
        }
        *ispsdefinite=psd;
        if (flag==DUAL_FACTOR){
          //DSDPEventLogEnd(ConeComputeS);
        } else if (flag==PRIMAL_FACTOR){
          //DSDPEventLogEnd(ConeComputeSS);
        }
      //}
      } // end of inlined function body of DSDPComputeSS
      if (psdefinite==DSDP_TRUE){
        dsdp->pstep=1.0;
        info=DSDPSaveYForX(dsdp,dsdp->mutarget,dsdp->pstep);DSDPCHKERR(info);
      } else {
        dsdp->pstep=0.0;
      }
  
      if (dsdp->usefixedrho==DSDP_TRUE){
        dsdp->rho=dsdp->rhon*dsdp->np;
        mutarget=(dsdp->ppobj-dsdp->ddobj)/(dsdp->rho);
        dsdp->pstep=0.5;
      } else {
        info = DSDPChooseBarrierParameter(dsdp,dsdp->mutarget,&dsdp->pstep,&mutarget);DSDPCHKERR(info);
        dsdp->rho=(dsdp->ppobj-dsdp->ddobj)/(mutarget);
      }
      DSDPEventLogEnd(dsdp->ptime);
      
      DSDPLogInfo(0,6,"Current mu=%4.8e, Target X with mu=%4.8e, Rho: %8.4e, Primal Step Length: %4.8f, pnorm: %4.8e\n",dsdp->mu,mutarget,dsdp->rho/dsdp->np,dsdp->pstep, dsdp->pnorm);
  
      /* Take Dual Step */
      /* Compute distance from chosen point on central path Pnorm */
      DSDPEventLogBegin(dsdp->dtime);
      info=DSDPComputeDY(dsdp,mutarget,dsdp->dy,&dsdp->pnorm); DSDPCHKERR(info);
      if (dsdp->pnorm<0.1){ mutarget/=10;  info=DSDPComputeDY(dsdp,mutarget,dsdp->dy,&dsdp->pnorm); DSDPCHKERR(info);}
  
      info=DSDPYStepLineSearch(dsdp, mutarget, 1.0, dsdp->dy);DSDPCHKERR(info);
      DSDPEventLogEnd(dsdp->dtime);
  
      maxattempts=dsdp->reuseM;
      if (dsdp->dstep<1 && dsdp->rgap<1e-5) maxattempts=0;
      if (dsdp->dstep<1e-13) maxattempts=0;
      if (dsdp->rgap<1e-6) maxattempts=0;
      if (maxattempts>dsdp->reuseM) maxattempts=dsdp->reuseM;
      for (attempt=0;attempt<maxattempts;attempt++){
        double cgtol=1e-6;
        if (attempt>0 && dsdp->pnorm < 0.1) break;
        if (attempt > 0 && dsdp->dstep<1e-4) break;
        if (dsdp->rflag) break;
        DSDPEventLogBegin(dsdp->ctime);
        DSDPLogInfo(0,2,"Reuse Matrix %d: Ddobj: %12.8e, Pnorm: %4.2f, Step: %4.2f\n",attempt,dsdp->ddobj,dsdp->pnorm,dsdp->dstep);
        info=DSDPInvertS(dsdp);DSDPCHKERR(info);
        info=DSDPComputeG(dsdp,dsdp->rhstemp,dsdp->rhs1,dsdp->rhs2);DSDPCHKERR(info);
        //Wei: inlined function body of DSDPComputeG
        {

        } // end of inlined function body DSDPComputeG
        if (dsdp->slestype==2 || dsdp->slestype==3){
  	if (dsdp->rflag){
	  //info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs1,dsdp->dy1,cgtol,&cg1);DSDPCHKERR(info);
	  {
          //int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
	    DSDPSchurMat MM=dsdp->M;
	    DSDPVec RHS = dsdp->rhs1;
	    DSDPVec X = dsdp->dy1;
            DSDPTruth *success = &cg1;

            int iter=0,n,info,maxit=10;
            double dd,ymax;
            DSDPCG *sles=dsdp->sles; 
            DSDPCGMat    CGM;
          
            info=DSDPEventLogBegin(dsdp->cgtime);
            info=DSDPVecZero(X);DSDPCHKERR(info);
            info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
            *success=DSDP_TRUE;
            if (0){
              maxit=0;
            } else if (dsdp->slestype==1){
          
              CGM.type=DSDPNoMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              cgtol*=1000;
              maxit=5;
          
            } else if (dsdp->slestype==2){
              CGM.type=DSDPUnfactoredMatrix;
              CGM.M=MM;
              CGM.Diag=sles->Diag;
              CGM.dsdp=dsdp;
              cgtol*=100;
              maxit=(int)sqrt(1.0*n)+10;
              if (maxit>20) maxit=20;
              info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
              /*
                info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
                info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
              */
              
            } else if (dsdp->slestype==3){
              CGM.type=DSDPFactoredMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              maxit=0;
              info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
              if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
              if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
                maxit=6;
              } else if (dsdp->rgap<1e-5){
                maxit=3;
              }
          
              //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
              //Wei: inlined function body of DSDPSchurMatSolve
              {
                
              //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                int info;
                //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                //Wei: inline function body DSDPSchurMatSolveM
                {  
                //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                //{
                  int info,n;
                  double *xx,*bb;
                  //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                  if (MM.dsdpops->matsolve){
                    info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                    info=DSDPVecZero(X);DSDPCHKERR(info);
                    info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                    info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                    info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                  } else {
                    //DSDPNoOperationError(MM);
          	  exit(-1);         // omit error handling
                  }
                  info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                  info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                  //info=DSDPEventLogEnd(hsolveevent);
                //}
                }   // end of inlined function body DSDPSchurMatSolveM
                info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
              //} // end of original function body
          
              } // end of inlined function DSDPSchurMatSolve
          
            } else if (dsdp->slestype==4){
              CGM.type=DSDPFactoredMatrix;
              CGM.M=MM;
              CGM.dsdp=dsdp;
              maxit=3;
              //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
              //Wei: inlined function body of DSDPSchurMatSolve
              {
                
              //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                int info;
                //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                //Wei: inline function body DSDPSchurMatSolveM
                {  
                //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                //{
                  int info,n;
                  double *xx,*bb;
                  //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                  if (MM.dsdpops->matsolve){
                    info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                    info=DSDPVecZero(X);DSDPCHKERR(info);
                    info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                    info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                    info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                    info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                  } else {
                    //DSDPNoOperationError(MM);
          	  exit(-1);         // omit error handling
                  }
                  info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                  info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                  //info=DSDPEventLogEnd(hsolveevent);
                //}
                }   // end of inlined function body DSDPSchurMatSolveM
                info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
              //} // end of original function body
          
              } // end of inlined function DSDPSchurMatSolve
            }
            if (n<maxit) maxit=n;
            
            info=DSDPConjugateGradient(CGM,X,RHS,
          			     sles->R,sles->BR,sles->P,sles->BP,
          			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);
          
            if (iter>=maxit) *success=DSDP_FALSE;
            info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
            if (dd<0) *success=DSDP_FALSE;
            info=DSDPEventLogEnd(dsdp->cgtime);
          //}


          } // end of inlined function body of DSDPCGSolve

	}
  	//info=DSDPCGSolve(dsdp,dsdp->M,dsdp->rhs2,dsdp->dy2,cgtol,&cg1);DSDPCHKERR(info);
	    {
            //int DSDPCGSolve(DSDP dsdp, DSDPSchurMat MM, DSDPVec RHS, DSDPVec X,double cgtol, DSDPTruth *success){
	      DSDPSchurMat MM=dsdp->M;
	      DSDPVec RHS = dsdp->rhs2;
	      DSDPVec X = dsdp->dy2;
              DSDPTruth *success = &cg1;

              int iter=0,n,info,maxit=10;
              double dd,ymax;
              DSDPCG *sles=dsdp->sles; 
              DSDPCGMat    CGM;
            
              info=DSDPEventLogBegin(dsdp->cgtime);
              info=DSDPVecZero(X);DSDPCHKERR(info);
              info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
              *success=DSDP_TRUE;
              if (0){
                maxit=0;
              } else if (dsdp->slestype==1){
            
                CGM.type=DSDPNoMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                cgtol*=1000;
                maxit=5;
            
              } else if (dsdp->slestype==2){
                CGM.type=DSDPUnfactoredMatrix;
                CGM.M=MM;
                CGM.Diag=sles->Diag;
                CGM.dsdp=dsdp;
                cgtol*=100;
                maxit=(int)sqrt(1.0*n)+10;
                if (maxit>20) maxit=20;
                info=DSDPVecSet(1.0,sles->Diag);DSDPCHKERR(info);
                /*
                  info=DSDPSchurMatGetDiagonal(MM,sles->Diag);DSDPCHKERR(info);
                  info=DSDPVecReciprocalSqrt(sles->Diag); DSDPCHKERR(info);
                */
                
              } else if (dsdp->slestype==3){
                CGM.type=DSDPFactoredMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                maxit=0;
                info=DSDPGetMaxYElement(dsdp,&ymax);DSDPCHKERR(info);
                if (ymax > 1e5 && dsdp->rgap<1e-1) maxit=3;
                if (0 && ymax > 1e5 && dsdp->rgap<1e-2){ 
                  maxit=6;
                } else if (dsdp->rgap<1e-5){
                  maxit=3;
                }
            
                //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
                //Wei: inlined function body of DSDPSchurMatSolve
                {
                  
                //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                  int info;
                  //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                  //Wei: inline function body DSDPSchurMatSolveM
                  {  
                  //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                  //{
                    int info,n;
                    double *xx,*bb;
                    //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                    if (MM.dsdpops->matsolve){
                      info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                      info=DSDPVecZero(X);DSDPCHKERR(info);
                      info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                      info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                      info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                    } else {
                      //DSDPNoOperationError(MM);
            	  exit(-1);         // omit error handling
                    }
                    info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                    info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                    //info=DSDPEventLogEnd(hsolveevent);
                  //}
                  }   // end of inlined function body DSDPSchurMatSolveM
                  info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                  info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
                //} // end of original function body
            
                } // end of inlined function DSDPSchurMatSolve
            
              } else if (dsdp->slestype==4){
                CGM.type=DSDPFactoredMatrix;
                CGM.M=MM;
                CGM.dsdp=dsdp;
                maxit=3;
                //info=DSDPSchurMatSolve(MM,RHS,X);DSDPCHKERR(info);
                //Wei: inlined function body of DSDPSchurMatSolve
                {
                  
                //int DSDPSchurMatSolve(DSDPSchurMat M, DSDPVec b, DSDPVec x){
                  int info;
                  //info=DSDPSchurMatSolveM(M,b,x);DSDPCHKERR(info);
                  //Wei: inline function body DSDPSchurMatSolveM
                  {  
                  //int DSDPSchurMatSolveM(DSDPSchurMat M, DSDPVec b, DSDPVec x)
                  //{
                    int info,n;
                    double *xx,*bb;
                    //info=DSDPEventLogBegin(hsolveevent);  //cannot track static
                    if (MM.dsdpops->matsolve){
                      info=DSDPVecGetArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecGetSize(X,&n); DSDPCHKERR(info);
                      info=DSDPVecZero(X);DSDPCHKERR(info);
                      info=DSDPVecGetArray(X,&xx); DSDPCHKERR(info);
                      info=(MM.dsdpops->matsolve)(MM.data,bb+1,xx+1,n-2); //DSDPChkMatError(MM,info);
                      info=DSDPVecRestoreArray(RHS,&bb); DSDPCHKERR(info);
                      info=DSDPVecRestoreArray(X,&xx); DSDPCHKERR(info);
                    } else {
                      //DSDPNoOperationError(MM);
            	  exit(-1);         // omit error handling
                    }
                    info=DSDPVecSetR(X,0.0);DSDPCHKERR(info);
                    info=DSDPVecSetC(X,0.0);DSDPCHKERR(info);
                    //info=DSDPEventLogEnd(hsolveevent);
                  //}
                  }   // end of inlined function body DSDPSchurMatSolveM
                  info=DSDPApplySMW(MM,RHS,X);DSDPCHKERR(info);
                  info=DSDPZeroFixedVariables(MM,X);DSDPCHKERR(info);
                //} // end of original function body
            
                } // end of inlined function DSDPSchurMatSolve
              }
              if (n<maxit) maxit=n;
              
              info=DSDPConjugateGradient(CGM,X,RHS,
            			     sles->R,sles->BR,sles->P,sles->BP,
            			     sles->TTT,cgtol,maxit,&iter);DSDPCHKERR(info);
            
              if (iter>=maxit) *success=DSDP_FALSE;
              info=DSDPVecDot(RHS,X,&dd);DSDPCHKERR(info);
              if (dd<0) *success=DSDP_FALSE;
              info=DSDPEventLogEnd(dsdp->cgtime);
            //}


            } // end of inlined function body of DSDPCGSolve
        }
        info=DSDPVecDot(dsdp->b,dsdp->dy1,&dd1);DSDPCHKERR(info);
        info=DSDPVecDot(dsdp->b,dsdp->dy2,&dd2);DSDPCHKERR(info);
        if (dd1>0 && dd2>0){
  	mutarget=DSDPMin(mutarget,dd1/dd2*dsdp->schurmu);
        }
        mutarget=mutarget*(dsdp->np/(dsdp->np+sqrt(dsdp->np)));	  
        info=DSDPComputeDY(dsdp,mutarget,dsdp->dy, &ppnorm);DSDPCHKERR(info);
        if (ppnorm<=0){ DSDPEventLogEnd(dsdp->ctime);  break; }
        dsdp->pnorm=ppnorm;
        //info=DSDPYStepLineSearch2(dsdp, mutarget, dsdp->dstep, dsdp->dy);DSDPCHKERR(info);
        //Wei: inlined function body of DSDPYStepLineSearch2
        {
        //int DSDPYStepLineSearch2(DSDP dsdp, double mutarget, double dstep0, DSDPVec dy){
	  double dstep0 = dsdp->dstep;    //passing param 3
          
          /* The merit function is the objective (DD) plus the barrier function */
          /* This line search is used in the corrector steps */
          int info, attempt, maxattempts=10;
          double dstep,newpotential,bdotdy,oldpotential,logdet;
          double maxmaxstep=0,steptol=1e-6;
          double a,b;
          DSDPTruth psdefinite;
          //info=DSDPComputeMaxStepLength(dsdp,dy,DUAL_FACTOR,&maxmaxstep);DSDPCHKERR(info);
          //Wei: inlined function body of DSDPComputeMaxStepLength 
          {
          //int DSDPComputeMaxStepLength(DSDP dsdp, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
            DSDPDualFactorMatrix flag = DUAL_FACTOR;           //passing param 3
            double* maxsteplength = &maxmaxstep;               //passing param 4
            int info,kk;
            double msteplength=1.0e30,conesteplength;
        
            if (flag==DUAL_FACTOR){
              //DSDPEventLogBegin(ConeMaxDStep);
            } else if (flag==PRIMAL_FACTOR){
              //DSDPEventLogBegin(ConeMaxPStep);
            }
            for (kk=0;kk<dsdp->ncones;kk++){
              DSDPEventLogBegin(dsdp->K[kk].coneid);
              conesteplength=1.0e20;
              //info=DSDPConeComputeMaxStepLength(dsdp->K[kk].cone,DY,flag,&conesteplength);DSDPCHKCONEERR(kk,info);
              //Wei:  inlined function body of DSDPConeComputeMaxStepLength
              {
              //int DSDPConeComputeMaxStepLength(DSDPCone K, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
                DSDPCone K=dsdp->K[kk].cone;       // passing param 1
                
                int info;
                double inner_conesteplength=1.0e20;
                conesteplength=1.0e30;
                if (K.dsdpops->conemaxsteplength){
                  info=K.dsdpops->conemaxsteplength(K.conedata,dsdp->dy,flag,&inner_conesteplength);//DSDPChkConeError(K,info);
                } else {
                  //DSDPNoOperationError(K);
          	exit(-1);                      //omit error handling
                }
                //*maxsteplength=conesteplength;
                conesteplength = inner_conesteplength;
              //}
          
              } // end of inlined function body of DSDPConeComputeMaxStepLength
              msteplength=DSDPMin(msteplength,conesteplength);
              DSDPEventLogEnd(dsdp->K[kk].coneid);
            }
            *maxsteplength=msteplength;
            if (flag==DUAL_FACTOR){
              //DSDPEventLogEnd(ConeMaxDStep);
            } else if (flag==PRIMAL_FACTOR){
              //DSDPEventLogEnd(ConeMaxPStep);
            }
          //} //original function body end
          }  // end of inlined function body of DSDPComputeMaxStepLength
          info=DSDPComputePotential2(dsdp,dsdp->y,mutarget, dsdp->logdet,&oldpotential);DSDPCHKERR(info);
          info=DSDPVecDot(dsdp->rhs,dsdp->dy,&bdotdy);DSDPCHKERR(info);
          dstep=DSDPMin(dstep0,0.95*maxmaxstep);
          if (dstep * dsdp->pnorm > dsdp->maxtrustradius) dstep=dsdp->maxtrustradius/dsdp->pnorm;
          DSDPLogInfo(0,8,"Full Dual StepLength %4.4e, %4.4e\n",maxmaxstep,dstep);
          for (psdefinite=DSDP_FALSE,attempt=0; attempt<maxattempts && psdefinite==DSDP_FALSE; attempt++){
            if (dstep < steptol) break;
            info=DSDPComputeNewY(dsdp,dstep,dsdp->ytemp);DSDPCHKERR(info);
            //info=DSDPComputeSS(dsdp,dsdp->ytemp,DUAL_FACTOR,&psdefinite);DSDPCHKERR(info);
            //Wei: inlined function body of DSDPComputeSS 
            {
            //int DSDPComputeSS(DSDP dsdp, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
              DSDPVec Y=dsdp->ytemp;
              DSDPDualFactorMatrix flag = DUAL_FACTOR;
              DSDPTruth *ispsdefinite = &psdefinite;

              int info,kk;
              DSDPTruth psd=DSDP_TRUE;
              if (flag==DUAL_FACTOR){
                //DSDPEventLogBegin(ConeComputeS);
              } else if (flag==PRIMAL_FACTOR){
                //DSDPEventLogBegin(ConeComputeSS);
              }
              for (kk=dsdp->ncones-1; kk>=0 && psd==DSDP_TRUE;kk--){
                DSDPEventLogBegin(dsdp->K[kk].coneid);
                //info=DSDPConeComputeS(dsdp->K[kk].cone,Y,flag,&psd); DSDPCHKCONEERR(kk,info);
                //Wei: inlined function body of DSDPConeComputeS
                {
                //int DSDPConeComputeS(DSDPCone K, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
                  DSDPCone K=dsdp->K[kk].cone;
                  DSDPTruth *ispsdefinite = &psd;
                  int info;
                  if (K.dsdpops->conecomputes){
                    info=K.dsdpops->conecomputes(K.conedata,Y,flag,ispsdefinite); //DSDPChkConeError(K,info);
                  } else {
                    //DSDPNoOperationError(K);
            	exit(-1);   //omit error handling
                  }
                //}
                } // end of inlined function body of DSDPConeComputeS
                DSDPEventLogEnd(dsdp->K[kk].coneid);
              }
              *ispsdefinite=psd;
              if (flag==DUAL_FACTOR){
                //DSDPEventLogEnd(ConeComputeS);
              } else if (flag==PRIMAL_FACTOR){
                //DSDPEventLogEnd(ConeComputeSS);
              }
            //}
            } // end of inlined function body of DSDPComputeSS
            
            if (psdefinite==DSDP_TRUE){
              info=DSDPComputeLogSDeterminant(dsdp,&logdet);DSDPCHKERR(info);
              info=DSDPComputePotential2(dsdp,dsdp->ytemp,mutarget,logdet,&newpotential);DSDPCHKERR(info);
              b=bdotdy; a=2*(newpotential-oldpotential+bdotdy*dstep)/(dstep*dstep);
              if (newpotential>oldpotential-0.1*dstep*bdotdy ){
        	DSDPLogInfo(0,2,"Not sufficient reduction. Reduce stepsize.  Step:: %4.4e\n",dstep);
        	psdefinite=DSDP_FALSE;
        	if (b/a<dstep && b/a>0){ dstep=b/a;} else { dstep=dstep/2; } 
              } 
            } else {
              dstep=dstep/2.0;
              DSDPLogInfo(0,2,"Dual Matrix not Positive Definite: Reduce step %4.4e",dstep);
            }
          } /* Hopefully, the initial step size works and only go through loop once */
          if (psdefinite==DSDP_TRUE && dstep>=steptol){
            info=DSDPSetY(dsdp,dstep,logdet,dsdp->ytemp);DSDPCHKERR(info);
          } else {
            info=DSDPSetY(dsdp,0,dsdp->logdet,dsdp->y);DSDPCHKERR(info);
          }
        //} //end of the original function body

        } // end of inlined function body of DSDPYStepLineSearch2
        DSDPEventLogEnd(dsdp->ctime);
      }
      if (attempt>0)dsdp->dstep=1.0;
      
      dsdp->mutarget=DSDPMin(dsdp->mu,mutarget);
      info=DSDPGetRR(dsdp,&dd1);DSDPCHKERR(info);
      if (dsdp->itnow==0 && dsdp->xmaker[0].pstep<1.0 && dd1> 0 && dsdp->pstep<1.0 && dsdp->goty0==DSDP_FALSE){
        info=DSDPResetY0(dsdp);DSDPCHKERR(info); continue;
        dsdp->goty0=DSDP_FALSE;
      }
      
    } /* End of Dual Scaling Algorithm */
    
  //} 

  } // end of inlined function body of DSDPSolveDynamicRho
  if (dsdp->pstep==1){info=DSDPRefineStepDirection(dsdp,dsdp->xmakerrhs,dsdp->xmaker[0].dy);DSDPCHKERR(info);}
  if (dsdp->pdfeasible==DSDP_PDUNKNOWN) dsdp->pdfeasible=DSDP_PDFEASIBLE;
  info=DSDPEventLogEnd(dsdp->solvetime);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPCallMonitors"
/*!
\fn int DSDPCallMonitors(DSDP dsdp,DMonitor dmonitor[], int ndmonitors);
\param dsdp solver
\param dmonitor array of monitors
\param ndmonitors number of monitors.
\brief Call the monitor routines.
 */
int DSDPCallMonitors(DSDP dsdp,DMonitor dmonitor[], int ndmonitors){
  int i,info;
  DSDPFunctionBegin;
  for (i=0; i<ndmonitors;i++){
    info=(dmonitor[i].monitor)(dsdp,dmonitor[i].monitorctx);  DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}
/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "DSDPCheckConvergence"
/*!
\fn int DSDPCheckConvergence(DSDP dsdp,DSDPTerminationReason *reason);
\param dsdp solver
\param reason termination reason
\brief Check for convergence and monitor solution
 */
int DSDPCheckConvergence(DSDP dsdp,DSDPTerminationReason *reason){
  int info;
  DSDPTruth unbounded;

  DSDPFunctionBegin;
  info = DSDPGetConicDimension(dsdp,&dsdp->np);DSDPCHKERR(info);
  dsdp->rgap=(dsdp->ppobj-dsdp->ddobj)/(1.0+fabs(dsdp->ppobj)+fabs(dsdp->ddobj));
  dsdp->pstepold=dsdp->pstep;
  if (dsdp->reason==CONTINUE_ITERATING){
    if (dsdp->itnow>0){
      info=DSDPCheckForUnboundedObjective(dsdp,&unbounded);DSDPCHKERR(info);
      if (unbounded==DSDP_TRUE){
	dsdp->pdfeasible=DSDP_UNBOUNDED;
	info=DSDPSetConvergenceFlag(dsdp,DSDP_CONVERGED); DSDPCHKERR(info); 
      }
    }
    if (dsdp->reason==CONTINUE_ITERATING){
      if (dsdp->muold<dsdp->mutarget && dsdp->pstep==1 && dsdp->dstep==1 && dsdp->rgap<1e-5){
	info=DSDPSetConvergenceFlag(dsdp,DSDP_NUMERICAL_ERROR); DSDPCHKERR(info);
	DSDPLogInfo(0,2,"DSDP Finished: Numerical issues: Increase in Barrier function. \n");}
      if (dsdp->itnow >= dsdp->maxiter){
	info=DSDPSetConvergenceFlag(dsdp,DSDP_MAX_IT); DSDPCHKERR(info);} 
      if (dsdp->Mshift>dsdp->maxschurshift){
	info = DSDPSetConvergenceFlag(dsdp,DSDP_INDEFINITE_SCHUR_MATRIX); DSDPCHKERR(info);
      }
    } 
    info=DSDPCallMonitors(dsdp,dsdp->dmonitor,dsdp->nmonitors);DSDPCHKERR(info);
    info=DSDPMonitorCones(dsdp,0); DSDPCHKERR(info);
  }
  dsdp->muold=dsdp->mutarget;
  info = DSDPStopReason(dsdp,reason); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}



/* ---------------------------------------------------------- */ 
#undef __FUNCT__  
#define __FUNCT__ "DSDPTakeDown"
/*!
\fn int DSDPTakeDown(DSDP dsdp);
\param dsdp solver
\brief Destroy internal data structures.
 */
int DSDPTakeDown(DSDP dsdp){

  int i,info;

  DSDPFunctionBegin;
  DSDPValid(dsdp);
  info = DSDPVecDestroy(&dsdp->rhs);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->rhs1);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->rhs2);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->rhstemp);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->y);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->ytemp);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->dy1);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->dy2);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->dy);DSDPCHKERR(info);
  for (i=0;i<MAX_XMAKERS;i++){
    info = DSDPVecDestroy(&dsdp->xmaker[i].y);DSDPCHKERR(info);
    info = DSDPVecDestroy(&dsdp->xmaker[i].dy);DSDPCHKERR(info);
    info = DSDPVecDestroy(&dsdp->xmaker[i].rhs);DSDPCHKERR(info);
  }
  info = DSDPVecDestroy(&dsdp->xmakerrhs);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->y0);DSDPCHKERR(info);
  info = DSDPVecDestroy(&dsdp->b);DSDPCHKERR(info);

  info = DSDPCGDestroy(&dsdp->sles);DSDPCHKERR(info);
  info = DSDPDestroyCones(dsdp);DSDPCHKERR(info);
  info = DSDPSchurMatDestroy(&dsdp->M);DSDPCHKERR(info);
  info = DSDPGetConicDimension(dsdp,&dsdp->np);DSDPCHKERR(info);
  dsdp->setupcalled=DSDP_FALSE;
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetDestroyRoutine(DSDP dsdp, int (*fd)(void*), void* ctx);

\brief Set a routine that will be called during DSDPDestroy().
\param dsdp the solver
\param fd function pointer
\param ctx pointer to structure.
\sa DSDPDestroy()
*/
int DSDPSetDestroyRoutine(DSDP dsdp, int (*fd)(void*), void* ctx){
  int nd=dsdp->ndroutines;
  if (nd<10){
    dsdp->droutine[nd].f=fd;
    dsdp->droutine[nd].ptr=ctx;
    dsdp->ndroutines++;
  } else {
    printf("TOO MANY Destroy routines\n");
    return 1;
  }
  return 0;
}


/*!
\fn int DSDPDestroy(DSDP dsdp)
\brief Free the internal data structures of the solver and
the cones associated with it.

\param dsdp the solver
\sa DSDPCreate()
\sa DSDPSolve()
\sa DSDPSetup()
\ingroup DSDPBasic
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPDestroy"
int DSDPDestroy(DSDP dsdp){
  int i,info;
  DSDPFunctionBegin;
  DSDPValid(dsdp);
  for (i=0;i<dsdp->ndroutines;i++){
    info=(*dsdp->droutine[i].f)(dsdp->droutine[i].ptr);DSDPCHKERR(info);
  }
  info=DSDPTakeDown(dsdp);DSDPCHKERR(info);
  DSDPFREE(&dsdp,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
