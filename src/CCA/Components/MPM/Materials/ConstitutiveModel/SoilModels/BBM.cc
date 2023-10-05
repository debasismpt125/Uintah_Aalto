/*

The MIT License

Copyright (c) 1997-2010 Center for the Simulation of Accidental Fires and
Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI),
University of Utah.

License for the specific language governing rights and limitations under
Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

*/


#include <CCA/Components/MPM/Materials/ConstitutiveModel/SoilModels/BBM.h>
#include <CCA/Components/MPM/Materials/ConstitutiveModel/SoilModels/BBMModel.h>
#include <CCA/Components/MPM/Materials/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/ProblemSpec/ProblemSpec.h>

//#include <Core/Containers/StaticArray.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/MinMax.h>

#include <sci_defs/uintah_defs.h>

#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <Core/Exceptions/InvalidValue.h>

//sets of external variables for the Sheng Mohr Coulomb algorithm by WTS. Some are redundant.
/*
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Those external variables are defined in the Mohr-Coulomb procedure. Therefore they are not defined again here
//if you remove those from there, you need to re-enable the lines below where the external variables are defined
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
//double ALFACHECK, ALFACHANGE, ALFARATIO, MAXITER, YIELDTOL, TOL_METHOD, INCREMENT_TYPE, BETA_FACT;
//double DRIFT_CORRECTION, EULER_ITERATIONS,CRITICAL_STEP_SIZE;
//double STEP_MAX, STEP_MIN, ERROR_DEF, USE_ERROR_STEP, MIN_DIVISION_SIZE;
//double INTEGRATION_TOL=0.0001;
//double LINES_SEC_NO=50;
//int MAX_LOOP, SOLUTION_ALGORITHM, ALGORITHM_TYPE;   //MAX_LOOP - number of steps to solve, SOLUTION_ALGORITHM - algorithm to use
//double USE_ERROR=0.5, SAVE_FOR_ERROR=1;    //these are values - 1st - when 'use' the additional error to speed up the calculations (here 30%)
// SAVEFORERROR says how greater accuracy the program should use; 1 means no change. 0.5 means 50% greater accuracy. (tol * 0.5) etc
//double CHGEPSINTOL=10e-9;
//double ADDTOLYIELD=0.8;
//double SUCTIONTOL=0.00000001;
//double TINY=1e-14;		//value used for example in checking whether we are not dividing by zero in CalcStressElast, volumetric strain must be also larger then tiny
//double PMIN=0.0001;			// value of minimum mean stress to calculate K in CalcStressElast
//int USE_NICE_SCHEME=0;

////////////////////////////////////////////////////////////////////////////////
// The following functions are found in fortran/*.F
//SUBROUTINE YENTA_CALC( NBLK, NINSV, DT, PROP,
  //   $                                   SIGARG, D, SVARG, USM   )

extern "C"{

#if defined( FORTRAN_UNDERSCORE_END )
#  define DMMCHK dmmchk_
#  define DIAMM_CALC diamm_calc_
#  define DMMRXV dmmrxv_
#elif defined( FORTRAN_UNDERSCORE_LINUX )
#  define DMMCHK dmmchk_
#  define DMMRXV dmmrxv_
#  define DIAMM_CALC diamm_calc__
#else // NONE
#  define DMMCHK dmmchk
#  define DIAMM_CALC diamm_calc
#  define  dmmrxv
#endif

//#define DMM_ANISOTROPIC
//#undef DMM_ANISOTROPIC

 //  void DMMCHK( double UI[], double UJ[], double UK[] );

/*
C
C***********************************************************************
C     REQUIRED MIG DATA CHECK ROUTINE
C     Checks validity of user inputs for DMM model.
C     Sets defaults for unspecified user input.
C     Adjusts user input to be self-consistent.
C
C     input
C     -----
C       UI: user input as read and stored by host code.
C
C       Upon entry, the UI array contains the user inputs EXACTLY
C       as read from the user.  These inputs must be ordered in the
C       UI array as indicated in the file kmmpnt.Blk.
C       See kmmpnt.Blk for parameter definitions and keywords.
C
C       DC: Not used with this model
C
C    Other output
C    ------------
C       GC: Not used with this model
C       DC: Not used with this model
C       Because GC and DC are not used, you may call this routine
C       with a line of the form "CALL DMMCHK(UI,UI,UI)"
*/

  // void DIAMM_CALC( int &nblk, int &ninsv, double &dt,
   //                                 double UI[], double stress[], double D[],
     //                               double svarg[], double &USM );
/*
C***********************************************************************
C
C     Description:
C           Drucker-Prager plasticity model with elastic strain induced
C           anisotropy.
C
C***********************************************************************
C
C     input arguments
C     ===============
C      NBLK       int                   Number of blocks to be processed
C      NINSV      int                   Number of internal state vars
C      DTARG      dp                    Current time increment
C      UI       dp,ar(nprop)            User inputs
C      D          dp,ar(6)              Strain increment
C
C     input output arguments
C     ======================
C      STRESS   dp,ar(6)                stress
C      SVARG    dp,ar(ninsv)            state variables
C
C     output arguments
C     ================
C      USM      dp                      uniaxial strain modulus
C
C***********************************************************************
C
C      stresss and strains, plastic strain tensors
C          11, 22, 33, 12, 23, 13
C
C***********************************************************************
*/

   void DMMRXV( double UI[], double UJ[], double UK[], int &nx, char* namea[],
                char* keya[], double rinit[], double rdim[], int iadvct[],
                int itype[] );
}
/*
C**********************************************************************
C     REQUESTED EXTRA VARIABLES FOR KAYENTA
C
C     This subroutine creates lists of the internal state variables
C     needed for DMM. This routine merely sends a
C     LIST of internal state variable requirements back to the host
C     code.   IT IS THE RESPONSIBILITY OF THE HOST CODE to loop over
C     the items in each list to actually establish necessary storage
C     and (if desired) set up plotting, restart, and advection
C     control for each internal state variable.
C
C     called by: host code after all input data have been checked
C
C     input
C     -----
C          UI = user input array
C          GC = unused for this model (placeholder)
C          DC = unused for this model (placeholder)
C
C     output
C     ------
C          NX = number of extra variables                    [DEFAULT=0]
C       NAMEA = single character array created from a string array
C               (called NAME) used locally in this routine to register
C               a descriptive name for each internal state variable.
C        KEYA = single character array created from a string array
C               (called KEY) used locally in this routine to register
C               a plot keyword for each internal state variable.
C          | Note: NAMEA and KEYA are created from the local variables |
C          | NAME and KEY by calls to the subroutine TOKENS, which     |
C          | is a SERVICE routine presumed to ALREADY exist within the |
C          | host code (we can provide this routine upon request).     |
C          | "NAME" is a fortran array of strings. "NAMEA" is a one    |
C          | dimensional array of single characters. For readability,  |
C          | most of this subroutine writes to the NAME array. Only at |
C          | the very end is NAME converted to NAMEA by calling a      |
C          | the utility routine TOKENS. The KEY array is similarly    |
C          | converted to KEYA.  These conversions are performed       |
C          | because host codes written in C or C++ are unable to      |
C          | process FORTRAN string arrays. Upon request, we can       |
C          | provide a utility routine that will convert BACK to       |
C          | FORTRAN string arrays if your host code is FORTRAN.       |
C          | Likewise, we can provide a C++ routine that will allow    |
C          | parsing the single-character arrays to convert them back  |
C          | to strings if your code is C++. Alternatively, you can    |
C          | simply ignore the NAMEA and KEYA outputs of this routine  |
C          | if your host code does not wish to establish plotting     |
C          | information.                                              |
C
C       RINIT = initial value for each ISV               [DEFAULT = 0.0]
C        RDIM = physical dimension exponents             [DEFAULT = 0.0]
C               This variable is dimensioned RDIM(7,*) for the 7 base
C               dimensions (and * for the number of extra variables):
C
C                      1 --- length
C                      2 --- mass
C                      3 --- time
C                      4 --- temperature
C                      5 --- discrete count
C                      6 --- electric current
C                      7 --- luminous intensity
C
C                Suppose, for example, that an ISV has units of stress.
C                Dimensionally, stress is length^(1) times mass^(-1)
C                times time^(-2). Therefore, this routine would return
C                1.0, -1.0, and -2.0 as the first three values of the
C                RDIM array. Host codes that work only in one unit
C                set (e.g., SI) typically ignore the RDIM output.
C
C      IADVCT = advection option                           [DEFAULT = 0]
C                    = 0 advect by mass-weighted average
C                    = 1 advect by volume-weighted average
C                    = 2 don't advect
C            The advection method will often be ignored by host codes.
C            It is used for Eulerian codes and for Lagrangian codes that
C            re-mesh (and therefore need guidance about how to "mix"
C            internal state variables). Note: a value of 2 implies that
C            the ISV is output only.
C
C        ITYPE = variable type (see migtionary preface)    [DEFAULT = 1]
C                  1=scalar
C                  6=2nd-order symmetric tensor
C        The component ordering for ITYPE=6 is 11, 22, 33, 12, 23, 31.
C        Consequently, the 11 component is the first one to be requested
C        in tensor lists, and its IFLAG is set to 6. To indicate that
C        subsequent ISVs are the remaining components of the same tensor,
C        the next five ISVs are given an IFLAG value of -6.
C        Host codes that don't change basis can ignore ITYPE.
*/
// End fortran functions.
////////////////////////////////////////////////////////////////////////////////
using std::cerr; using namespace Uintah;

BBM::BBM(ProblemSpecP& ps,MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  d_NBASICINPUTS=34;
  d_NMGDC=13;

// Total number of properties
  d_NDMMPROP=d_NBASICINPUTS+d_NMGDC;

  // pre-initialize all of the user inputs to zero.
  for(int i = 0; i<d_NDMMPROP; i++){
     UI[i] = 0.;
  }
  // Read model parameters from the input file
  getInputParameters(ps);

  // Check that model parameters are valid and allow model to change if needed


  //DMMCHK(UI,UI,&UI[d_NBASICINPUTS]);
  CheckModel (UI);  //note no check implemented for BBM

  //Create VarLabels for GeoModel internal state variables (ISVs)
  int nx;


 // DMMRXV( UI, UI, UI, nx, namea, keya, rinit, rdim, iadvct, itype );

  nx=27;
     for (int i=0;i<nx;i++)
    {
        rinit[i]=UI[i];
        //cerr<<" UI["<<i<<"]="<<UI[i];
    }

  d_NINSV=nx;
    //cerr << "d_NINSV = " << d_NINSV <<endl;
    //getchar();
  initializeLocalMPMLabels();
}
#if 0
BBM::BBM(const BBM* cm) : ConstitutiveModel(cm)
{
  for(int i=0;i<d_NDMMPROP;i++){
    UI[i] = cm->UI[i];
  }

  //Create VarLabels for Diamm internal state variables (ISVs)
  initializeLocalMPMLabels();
}
#endif
BBM::~BBM()
{
   for (unsigned int i = 0; i< ISVLabels.size();i++){
     VarLabel::destroy(ISVLabels[i]);
   }
}

void BBM::outputProblemSpec(ProblemSpecP& ps,bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type","BBM");
  }


/*
     G - shear modulus,
     Kappa - swelling index for changes in mean net stress
     KappaS - Swelling index for changes in suction
     p_atm - Atmospheric pressure
     Lambda -Slope of Normal Compression Line at zero suction
     r - parameter controlling ration of NCL slopes at s->inf and s=0
     BetaBBM - parameter controlling variation of NCL slope with suction
     p_c - reference pressure
     k - parameter controlling cohesion increase with suction
     M - slope of critical state line
     N_zero - specific volume at reference pressure
     Non_Associated_Flow - 0-no, 1 yes
    ------- hardening paramters / variables
    p_star - hardening parameter
    s - suction
*/
  //PMIN=0.0; //need to set initial stress to a relatively high value
  //'cout<<"WARNING: due to requirement of stress being positive, initial stress of particles has been set to 1.0. This means that the hardening value etc set in BBM must be >>1 for the algorithm to work correctly."<<endl;


  cm_ps->appendElement("G_BBM",UI[0]);   //  shear modulus (stress)
  cm_ps->appendElement("Kappa",UI[1]);
  cm_ps->appendElement("KappaS",UI[2]);
  cm_ps->appendElement("p_atm",UI[3]);
  cm_ps->appendElement("Lambda",UI[4]);
  cm_ps->appendElement("r",UI[5]);
  cm_ps->appendElement("BetaBBM",UI[6]);
  cm_ps->appendElement("p_c",UI[7]);
  cm_ps->appendElement("k",UI[8]);
  cm_ps->appendElement("M",UI[9]);
  cm_ps->appendElement("Non_Associated_Flow",UI[10]);
  cm_ps->appendElement("N_zero",UI[11]);
  cm_ps->appendElement("p_star",UI[12]);
  cm_ps->appendElement("suction",UI[13]);  //
  cm_ps->appendElement("p_initial",UI[14]);  //
  cm_ps->appendElement("BulkModulusK",UI[15]);  //
  cm_ps->appendElement("SpecificVolume",UI[16]);  //

  // below: unused for the model at the moment; left as they were for the DIAMM model
/*
   double PStar=UI[12],NZero=UI[11],KappaS=UI[2],suction=UI[13],pc=UI[7];
   double PAtmos=UI[3],KappaP=UI[1],Beta=UI[6],r=UI[5],LambdaZero=UI[4];

   double   Ns=NZero-KappaS*log((suction+PAtmos)/PAtmos);
   double   LambdaS=LambdaZero*((1-r)*exp(-Beta*suction)+r);
   double  PZero=(LambdaZero-KappaP)/(LambdaS-KappaP);
	  PZero=pc*pow((PStar/pc),PZero);
      double nu0=Ns-LambdaS*log(PZero/pc); //assuming that on the yield surface
      double K=PZero*nu0/KappaP; //assuming on the yield surface, updated later
      //note: with initial p it can be better computed
    UI[15]=K;
    UI[16]=nu0;

*/

  cm_ps->appendElement("PMIN",UI[17]);  //
  cm_ps->appendElement("Accuracy",UI[18]);  //
  cm_ps->appendElement("CriticalSubstepSize",UI[19]);  //
  cm_ps->appendElement("DriftCorrection",UI[20]);  //
  cm_ps->appendElement("IntegrationAlgorithm",UI[21]);  //

  cm_ps->appendElement("UseWaterRetention",UI[22]);  //
  cm_ps->appendElement("WR_Param1",UI[23]);  //
  cm_ps->appendElement("WR_Param2",UI[24]);  //
  cm_ps->appendElement("WR_Param3",UI[25]);  //
  cm_ps->appendElement("WR_Param4",UI[26]);  //


  cm_ps->appendElement("IDK",UI[27]);//
  cm_ps->appendElement("IDG",UI[28]);//

  cm_ps->appendElement("A4PF",UI[29]);//

  cm_ps->appendElement("TQC",UI[30]);//
  cm_ps->appendElement("F1",UI[31]);//

  cm_ps->appendElement("TEST",UI[32]);//
  cm_ps->appendElement("DEJAVU",UI[33]);//

  int dcprop=d_NBASICINPUTS-1;
  cm_ps->appendElement("DC1",UI[dcprop+1]);//
  cm_ps->appendElement("DC2",UI[dcprop+2]);//
  cm_ps->appendElement("DC3",UI[dcprop+3]);//
  cm_ps->appendElement("DC4",UI[dcprop+4]);//
  cm_ps->appendElement("DC5",UI[dcprop+5]);//
  cm_ps->appendElement("DC6",UI[dcprop+6]);//
  cm_ps->appendElement("DC7",UI[dcprop+7]);//
  cm_ps->appendElement("DC8",UI[dcprop+8]);//
  cm_ps->appendElement("DC9",UI[dcprop+9]);//
  cm_ps->appendElement("DC10",UI[dcprop+10]);//
  cm_ps->appendElement("DC11",UI[dcprop+11]);//
  cm_ps->appendElement("DC12",UI[dcprop+12]);//
  cm_ps->appendElement("DC13",UI[dcprop+13]);//

/*
  cerr<<"outputProblemSpec procedure: specification:";
  for (int i=0; i<17; i++) cerr<<(UI[i])<<", ";
  cerr<<endl;
  */
}

BBM* BBM::clone()
{
  return scinew BBM(*this);
}

void BBM::initializeCMData(const Patch* patch,
                               const MPMMaterial* matl,
                               DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
    std::vector<ParticleVariable<double> > ISVs(d_NINSV+1);
  cout << "In initializeCMData" << endl;
  for(int i=0;i<d_NINSV;i++){
    new_dw->allocateAndPut(ISVs[i],ISVLabels[i], pset);
    ParticleSubset::iterator iter = pset->begin();
    for(;iter != pset->end(); iter++){
      ISVs[i][*iter] = rinit[i];
    }
  }

  computeStableTimestep(patch, matl, new_dw);
}


void BBM::addParticleState(std::vector<const VarLabel*>& from,
                               std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  for(int i=0;i<d_NINSV;i++){
    from.push_back(ISVLabels[i]);
    to.push_back(ISVLabels_preReloc[i]);
  }
}

void BBM::computeStableTimestep(const Patch* patch,
                                    const MPMMaterial* matl,
                                    DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int dwi = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  constParticleVariable<double> pmass, pvolume;
  constParticleVariable<Vector> pvelocity;
  //constParticleVariable<Matrix3> pstress,deformationGradient;

  new_dw->get(pmass,     lb->pMassLabel,     pset);
  new_dw->get(pvolume,   lb->pVolumeLabel,   pset);
  new_dw->get(pvelocity, lb->pVelocityLabel, pset);
  //new_dw->get(pstress, lb->pStressLabel, pset);
  //new_dw->get(deformationGradient, lb->pDeformationMeasureLabel, pset);
  //Matrix3 tensorR, tensorU,tensorSig;
  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12,1.e-12,1.e-12);
   std::vector<constParticleVariable<double> > ISVs(d_NINSV+1);
    for(int i=0;i<d_NINSV;i++){
      new_dw->get(ISVs[i],           ISVLabels[i],                 pset);
    }
  //
  //double bulk, meanstress,nu0,Ns,PZero,suction,PStar;
  //double G = UI[0], KappaP=UI[1],KappaS=UI[2],PAtmos=UI[3]; //modified: K=UI[1], G=UI[0]
  //double LambdaZero=UI[4],r=UI[5],Beta=UI[6],pc=UI[7],NZero=UI[11];
  double bulk,G;
  //double svarg[d_NINSV];

  for(ParticleSubset::iterator iter = pset->begin();iter != pset->end();iter++){
     particleIndex idx = *iter;
      bulk=ISVs[15][idx];
      G=ISVs[0][idx];
     // Compute wave speed at each particle, store the maximum
     c_dil = sqrt((bulk + 4.*G/3.)*pvolume[idx]/pmass[idx]);
     WaveSpeed=Vector(Max(c_dil+fabs(pvelocity[idx].x()),WaveSpeed.x()),
                      Max(c_dil+fabs(pvelocity[idx].y()),WaveSpeed.y()),
                      Max(c_dil+fabs(pvelocity[idx].z()),WaveSpeed.z()));
  }
  //UI[14]=matl->getInitialDensity();
  //UI[15]=matl->getRoomTemperature();
  //UI[14]=bulk/matl->getInitialDensity();  ??tim
  //UI[19]=matl->getInitialCv();
  WaveSpeed = dx/WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}



void BBM::computeStressTensor(const PatchSubset* patches,
                                  const MPMMaterial* matl,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
    //cout<<"Entered computeStressTensor Procedure"<<endl;
  double rho_orig = matl->getInitialDensity();
  for(int p=0;p<patches->size();p++){
    double se = 0.0;
    const Patch* patch = patches->get(p);

    //ParticleInterpolator* interpolator = flag->d_interpolator->clone(patch);
    //vector<IntVector> ni(interpolator->size());
    //vector<Vector> d_S(interpolator->size());
    //vector<double> S(interpolator->size());

    Matrix3 Identity,zero(0.),One(1.);Identity.Identity();;
    double c_dil=0.0;
    Vector WaveSpeed(1.e-12,1.e-12,1.e-12);
    Vector dx = patch->dCell();
    //double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    int dwi = matl->getDWIndex();
    // Create array for the particle position
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    constParticleVariable<Point> px, pxnew;
    constParticleVariable<Matrix3> deformationGradient, pstress;
    ParticleVariable<Matrix3> pstress_new;
    constParticleVariable<Matrix3> deformationGradient_new, velGrad;
    constParticleVariable<double> pmass, pvolume, ptemperature;
    constParticleVariable<double> pvolume_new;
    constParticleVariable<Vector> pvelocity;
    constParticleVariable<Matrix3> psize;

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

 //   Ghost::GhostType  gac   = Ghost::AroundCells;

    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    old_dw->get(px,                  lb->pXLabel,                  pset);
    old_dw->get(pstress,             lb->pStressLabel,             pset);
    old_dw->get(pmass,               lb->pMassLabel,               pset);
    old_dw->get(pvolume,             lb->pVolumeLabel,             pset);
    old_dw->get(pvelocity,           lb->pVelocityLabel,           pset);
    old_dw->get(ptemperature,        lb->pTemperatureLabel,        pset);
    old_dw->get(deformationGradient, lb->pDeformationMeasureLabel, pset);
	new_dw->get(pvolume_new,		 lb->pVolumeLabel_preReloc,    pset);
	new_dw->get(pxnew,				 lb->pXLabel_preReloc,		   pset);


    std::vector<constParticleVariable<double> > ISVs(d_NINSV+1);
    for(int i=0;i<d_NINSV;i++){
      old_dw->get(ISVs[i],           ISVLabels[i],                 pset);
    }

   // new_dw->get(gvelocity,lb->gVelocityStarLabel, dwi,patch, gac, NGN);

    ParticleVariable<double> pdTdt,p_q;

    new_dw->allocateAndPut(pstress_new,     lb->pStressLabel_preReloc,   pset);
    new_dw->allocateAndPut(pdTdt,           lb->pdTdtLabel,              pset);
    new_dw->allocateAndPut(p_q,             lb->p_qLabel_preReloc,       pset);
    new_dw->get(deformationGradient_new,
                                 lb->pDeformationMeasureLabel_preReloc,  pset);
    new_dw->get(velGrad,         lb->pVelGradLabel_preReloc,             pset);

    std::vector<ParticleVariable<double> > ISVs_new(d_NINSV+1);
    for(int i=0;i<d_NINSV;i++){
      new_dw->allocateAndPut(ISVs_new[i],ISVLabels_preReloc[i], pset);
    }

    for(ParticleSubset::iterator iter = pset->begin();
                                        iter != pset->end(); iter++){
      particleIndex idx = *iter;

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;
      // Initialize velocity gradient
      //velGrad.set(0.0);

      /*if(!flag->d_axisymmetric){
        // Get the node indices that surround the cell
        interpolator->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],deformationGradient[idx]);

        computeVelocityGradient(velGrad,ni,d_S,oodx,gvelocity);

      } else {  // axi-symmetric kinematics
        // Get the node indices that surround the cell
        interpolator->findCellAndWeightsAndShapeDerivatives(px[idx],ni,S,d_S,
							    psize[idx],deformationGradient[idx]);
        // x -> r, y -> z, z -> theta
        computeAxiSymVelocityGradient(velGrad,ni,d_S,S,oodx,gvelocity,px[idx]);
      }*/

      // Calculate rate of deformation D, and deviatoric rate DPrime,
     Matrix3 D = (velGrad[idx] + velGrad[idx].Transpose())*.5;
    //Matrix3 D = (velGrad + velGrad.Transpose())*.5;

      double J = deformationGradient[idx].Determinant();
      // Check 1: Look at Jacobian
      if (!(J > 0.0)) {
        cerr << getpid() ;
        constParticleVariable<long64> pParticleID;
        old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
        cerr << "**ERROR** Negative Jacobian of deformation gradient"
             << " in particle " << pParticleID[idx] << endl;
        cerr << "l = " << velGrad[idx] << endl;
        cerr << "F_old = " << deformationGradient[idx] << endl;
        cerr << "F_new = " << deformationGradient_new[idx] << endl;
        cerr << "J = " << J << endl;
        throw InternalError("Negative Jacobian",__FILE__,__LINE__);
      }
    //  pvolume_new[idx]=Jinc*pvolume[idx];

      // Compute the local sound speed
      double rho_cur = rho_orig/J;

      // NEED TO FIND R
      Matrix3 tensorR, tensorU;

      // Look into using Rebecca's PD algorithm
      deformationGradient_new[idx].polarDecompositionRMB(tensorU, tensorR);

      // This is the previous timestep Cauchy stress
      // unrotated tensorSig=R^T*pstress*R
      Matrix3 tensorSig = (tensorR.Transpose())*(pstress[idx]*tensorR);

      // Load into 1-D array for the fortran code
      double sigarg[6];
      sigarg[0]=tensorSig(0,0);
      sigarg[1]=tensorSig(1,1);
      sigarg[2]=tensorSig(2,2);
      sigarg[3]=tensorSig(0,1);
      sigarg[4]=tensorSig(1,2);
      sigarg[5]=tensorSig(2,0);

      // UNROTATE D: S=R^T*D*R
      D=(tensorR.Transpose())*(D*tensorR);

      // Load into 1-D array for the fortran code
      double Darray[6];
      Darray[0]=D(0,0);
      Darray[1]=D(1,1);
      Darray[2]=D(2,2);
      Darray[3]=D(0,1);
      Darray[4]=D(1,2);
      Darray[5]=D(2,0);
      double svarg[d_NINSV];
      double USM=9e99;
      double dt = delT;
      int nblk = 1;

      // Load ISVs into a 1D array for fortran code
      for(int i=0;i<d_NINSV;i++){
        svarg[i]=ISVs[i][idx];
        }

// Calling the external model here
      //DIAMM_CALC(nblk, d_NINSV, dt, UI, sigarg, Darray, svarg, USM);
      CalculateStress (nblk, d_NINSV, dt, UI, sigarg, Darray, svarg, USM);


//svarg contain new internal state variables for given particle



      // Unload ISVs from 1D array into ISVs_new
      for(int i=0;i<d_NINSV;i++){
        ISVs_new[i][idx]=svarg[i];
      }

      // This is the Cauchy stress, still unrotated
     //additional (and redundant) check


      for (int i=0; i<6; i++) if (!isfinite(sigarg[i]))
    {
      cerr<<endl<<endl<<"Output stress ["<<i<<"] is not finite. Set to 0.001"<<endl<<endl<<endl;
      sigarg[i]=-0.001;
    }

      tensorSig(0,0) = sigarg[0];
      tensorSig(1,1) = sigarg[1];
      tensorSig(2,2) = sigarg[2];
      tensorSig(0,1) = sigarg[3];
      tensorSig(1,0) = sigarg[3];
      tensorSig(2,1) = sigarg[4];
      tensorSig(1,2) = sigarg[4];
      tensorSig(2,0) = sigarg[5];
      tensorSig(0,2) = sigarg[5];


      // ROTATE pstress_new: S=R*tensorSig*R^T
      pstress_new[idx] = (tensorR*tensorSig)*(tensorR.Transpose());
      //cout << pstress_new[idx] << endl;

#if 0
      cout << pstress_new[idx] << endl;
#endif

      c_dil = sqrt(USM/rho_cur);

      // Compute the strain energy for all the particles
      Matrix3 AvgStress = (pstress_new[idx] + pstress[idx])*.5;

      double e = (D(0,0)*AvgStress(0,0) +
                  D(1,1)*AvgStress(1,1) +
                  D(2,2)*AvgStress(2,2) +
              2.*(D(0,1)*AvgStress(0,1) +
                  D(0,2)*AvgStress(0,2) +
                  D(1,2)*AvgStress(1,2))) * pvolume_new[idx]*delT;

      se += e;

      // Compute wave speed at each particle, store the maximum
      Vector pvelocity_idx = pvelocity[idx];
      WaveSpeed=Vector(Max(c_dil+fabs(pvelocity_idx.x()),WaveSpeed.x()),
                       Max(c_dil+fabs(pvelocity_idx.y()),WaveSpeed.y()),
                       Max(c_dil+fabs(pvelocity_idx.z()),WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificial_viscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z())/3.0;
        double c_bulk = sqrt(UI[1]/rho_cur);
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    }  // end loop over particles


    WaveSpeed = dx/WaveSpeed;
    double delT_new = WaveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se),     lb->StrainEnergyLabel);
    }

    //delete interpolator;
  }
}

void BBM::carryForward(const PatchSubset* patches,
                           const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    int dwi = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    std::vector<constParticleVariable<double> > ISVs(d_NINSV+1);
    std::vector<ParticleVariable<double> > ISVs_new(d_NINSV+1);

    for(int i=0;i<d_NINSV;i++){
      old_dw->get(ISVs[i],ISVLabels[i], pset);
      new_dw->allocateAndPut(ISVs_new[i],ISVLabels_preReloc[i], pset);
      ISVs_new[i].copyData(ISVs[i]);
  }

    // Don't affect the strain energy or timestep size
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.),     lb->StrainEnergyLabel);
    }
  }

}

void BBM::addInitialComputesAndRequires(Task* task,
                                            const MPMMaterial* matl,
                                            const PatchSet* ) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();

  cout << "In add InitialComputesAnd" << endl;

  // Other constitutive model and input dependent computes and requires
  for(int i=0;i<d_NINSV;i++){
    task->computes(ISVLabels[i], matlset);
  }
}

void BBM::addComputesAndRequires(Task* task,
                                     const MPMMaterial* matl,
                                     const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);

  // Computes and requires for internal state data
  for(int i=0;i<d_NINSV;i++){
    task->requires(Task::OldDW, ISVLabels[i],          matlset, Ghost::None);
    task->computes(             ISVLabels_preReloc[i], matlset);
  }
}

void BBM::addComputesAndRequires(Task*,
                                     const MPMMaterial*,
                                     const PatchSet*,
                                     const bool ) const
{
}

double BBM::computeRhoMicroCM(double pressure,
                                  const double p_ref,
                                  const MPMMaterial* matl,
                                  double temperature,
                                  double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double p_gauge = pressure - p_ref;
  double rho_cur;
  double bulk = UI[15];  //note that it is *very* inaccurate, but otherwise lots of computations
  //Also, that procedure seems to be unused - may only be called from somewhere else in the code
  //So, maybe this is not a big deal... that we have the density badly [linearly) updated...

  rho_cur = rho_orig/(1-p_gauge/bulk);

  return rho_cur;

#if 1
  cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR BBM" << endl;
#endif
}

void BBM::computePressEOSCM(double rho_cur, double& pressure,
                                double p_ref,
                                double& dp_drho,      double& tmp,
                                const MPMMaterial* matl,
                                double temperature)
{

  double bulk = UI[15];

  //note that it is *very* inaccurate, but otherwise lots of computations
  //Also, that procedure seems to be unused - may only be called from somewhere else in the code
  //So, maybe this is not a big deal... that we have the density badly [linearly) updated...

  double rho_orig = matl->getInitialDensity();

  double p_g = bulk*(1.0 - rho_orig/rho_cur);
  pressure = p_ref + p_g;
  dp_drho  = bulk*rho_orig/(rho_cur*rho_cur);
  tmp = bulk/rho_cur;  // speed of sound squared

#if 1
  cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR BBM" << endl;
#endif
}

double BBM::getCompressibility()
{
return 1.0/UI[15];

//return 1; //experimental:1 // found to be irrelevant
}

void
BBM::getInputParameters(ProblemSpecP& ps)
{
  ps->getWithDefault("G_BBM",UI[0],0.0);              // shear modulus
  ps->getWithDefault("Kappa",UI[1],0.0);
  ps->getWithDefault("KappaS",UI[2],0.0);
  ps->getWithDefault("p_atm",UI[3],0.0);
  ps->getWithDefault("Lambda",UI[4],0.0);
  ps->getWithDefault("r",UI[5],0.0);
  ps->getWithDefault("BetaBBM",UI[6],0.0);
  ps->getWithDefault("p_c",UI[7],0.0);
  ps->getWithDefault("k",UI[8],0.0);
  ps->getWithDefault("M",UI[9],0.0);
  ps->getWithDefault("Non_Associated_Flow",UI[10],0);
  ps->getWithDefault("N_zero",UI[11],0.0);
  ps->getWithDefault("p_star",UI[12],0.0);
  ps->getWithDefault("suction",UI[13],0.0);  //
  ps->getWithDefault("p_initial",UI[14],0.0);  //
  ps->getWithDefault("BulkModulusK",UI[15],0.0);  //
  ps->getWithDefault("SpecificVolume",UI[16],0.0);  //
  ps->getWithDefault("PMIN",UI[17],0.0);  //
  ps->getWithDefault("Accuracy",UI[18],0.0);  //
  ps->getWithDefault("CriticalSubstepSize",UI[19],0.0);  //
  ps->getWithDefault("DriftCorrection",UI[20],0.0);  //
  ps->getWithDefault("IntegrationAlgorithm",UI[21],0.0);  //
  ps->getWithDefault("UseWaterRetention",UI[22],0.0);  //
  ps->getWithDefault("WR_Param1",UI[23],0.0);  //
  ps->getWithDefault("WR_Param2",UI[24],0.0);  //
  ps->getWithDefault("WR_Param3",UI[25],0.0);//
  ps->getWithDefault("WR_Param4",UI[26],0.0);//




  ps->getWithDefault("IDK",UI[27],0.0);//
  ps->getWithDefault("IDG",UI[28],0.0);//
  ps->getWithDefault("A2PF",UI[29],0.0);//
  ps->getWithDefault("TQC",UI[30],0.0);//
  ps->getWithDefault("F1",UI[31],0.0);//
  ps->getWithDefault("TEST",UI[32],0.0);//
  ps->getWithDefault("DEJAVU",UI[33],0.0);//

  ps->getWithDefault("DC1",UI[34],0.0);//
  ps->getWithDefault("DC2",UI[35],0.0);//
  ps->getWithDefault("DC3",UI[36],0.0);//
  ps->getWithDefault("DC4",UI[37],0.0);//
  ps->getWithDefault("DC5",UI[38],0.0);//
  ps->getWithDefault("DC6",UI[39],0.0);//
  ps->getWithDefault("DC7",UI[40],0.0);//
  ps->getWithDefault("DC8",UI[41],0.0);//
  ps->getWithDefault("DC9",UI[42],0.0);//
  ps->getWithDefault("DC10",UI[43],0.0);//
  ps->getWithDefault("DC11",UI[44],0.0);//
  ps->getWithDefault("DC12",UI[45],0.0);//
  ps->getWithDefault("DC13",UI[46],0.0);//
}

void
BBM::initializeLocalMPMLabels()
{
  vector<string> ISVNames;

  ISVNames.push_back("G_BBM+");
  ISVNames.push_back("Kappa+");
  ISVNames.push_back("KappaS+");
  ISVNames.push_back("p_atm");
  ISVNames.push_back("Lambda+");
  ISVNames.push_back("r+");
  ISVNames.push_back("BetaBBM+");
  ISVNames.push_back("p_c+");
  ISVNames.push_back("k+");
  ISVNames.push_back("M+");
  ISVNames.push_back("Naf");
  ISVNames.push_back("NZero+");
  ISVNames.push_back("P_Star");
  ISVNames.push_back("suction");
  ISVNames.push_back("pinitial");
  ISVNames.push_back("BulkModulusT");
  ISVNames.push_back("SpecificVol");
  ISVNames.push_back("MeanStressP");
  ISVNames.push_back("DeviatorStressQ");
  ISVNames.push_back("DegreeOfSaturation");
  ISVNames.push_back("EXX");
  ISVNames.push_back("EYY");
  ISVNames.push_back("EZZ");
  ISVNames.push_back("EXY");
  ISVNames.push_back("EYZ");
  ISVNames.push_back("EXZ");
  ISVNames.push_back("EJ2");
 // ISVNames.push_back("PSTAR");


  for(int i=0;i<d_NINSV;i++){
    ISVLabels.push_back(VarLabel::create(ISVNames[i],
                          ParticleVariable<double>::getTypeDescription()));
    ISVLabels_preReloc.push_back(VarLabel::create(ISVNames[i]+"+",
                          ParticleVariable<double>::getTypeDescription()));
  }

}

//CODE ADDED BY WTS FOR SHENGMOHRCOULOMB BELOW

void BBM::CalculateStress (int &nblk, int &ninsv, double &dt,
                                    double UI[], double stress[], double D[],
                                    double svarg[], double &USM)


/*
copied from DIAMM, giving the input required
The procedure calculates stress for the Mohr-Coulomb model, with several
additional options. Those include the variation of the flow rule,
several choices of Mohr-Coulomb like yield surfaces and the dependency
on the strain rate. Most of the above-described features are 'work in
progress'

The list of arguments below is from the original DIAMM model. The arguments
are a bit adjusted to the Mohr-Coulomb, though in spirit they remain same.

int &nblk, int &ninsv, double &dt,double UI[], double stress[], double D[],
double svarg[], double &USM
C
C     input arguments
C     ===============
C      NBLK       int                   Number of blocks to be processed
C      NINSV      int                   Number of internal state vars
C      DTARG      dp                    Current time increment
C      UI       dp,ar(nprop)            User inputs
C      D          dp,ar(6)              Strain increment
C
C     input output arguments
C     ======================
C      STRESS   dp,ar(6)                stress
C      SVARG    dp,ar(ninsv)            state variables
C
C     output arguments
C     ================
C      USM      dp                      uniaxial strain modulus
C
C***********************************************************************
C
C      stresss and strains, plastic strain tensors
C          11, 22, 33, 12, 23, 13
C
C***********************************************************************

*/

{

//this is slightly slow as each model used needs to be declared, but on the other hand allows for keeping things clean

BBMModel Model;


if (nblk!=1) cerr<<"Mohr-Coulomb model may only be used with nblk equal to 1. Results obtained are incorrect."<<endl;

/*
ps->getWithDefault("G",UI[0],0.0);              // shear modulus
  ps->getWithDefault("Kappa",UI[1],0.0);
  ps->getWithDefault("KappaS",UI[2],0.0);
  ps->getWithDefault("p_atm",UI[3],0.0);
  ps->getWithDefault("Lambda",UI[4],0.0);
  ps->getWithDefault("r",UI[5],0.0);
  ps->getWithDefault("BetaBBM",UI[6],0.0);
  ps->getWithDefault("p_c",UI[7],0.0);
  ps->getWithDefault("k",UI[8],0.0);
  ps->getWithDefault("M",UI[9],0.0);
  ps->getWithDefault("Non_Associated_Flow",UI[10],0);
  ps->getWithDefault("N_zero",UI[11],0.0);
  ps->getWithDefault("p_star",UI[12],0.0);
  ps->getWithDefault("suction",UI[13],0.0);  //
  ps->getWithDefault("QuickNu",0.0);  //
  ps->getWithDefault("QuickK",0.0);  //
*/


for (int i=0; i<23; i++)
{
    if (!isfinite(svarg[i]))
    {
      cerr<<endl<<endl<<"svarg ["<<i<<"] is not finite"<<endl<<endl<<endl;
      getchar();
    }
}


double G=svarg[0]; //shear modulus [stress units]
double KappaP=svarg[1];
double KappaS=svarg[2];
double PAtmos=svarg[3];
double Lambda=svarg[4];
double r=svarg[5];
double BetaBBM=svarg[6];
double pc=svarg[7];
double k=svarg[8];
double M=svarg[9];
bool NonAssociated=svarg[10];
double NZero=svarg[11];
double PStar=svarg[12];
double suction=svarg[13];
double initial_p=svarg[14];
double QuickNu=svarg[16];
double QuickK=svarg[15];
double UseWaterRetention=svarg [23];

//Water Retention Parameters
double WTRParam[5];
WTRParam[0]=svarg[19];
WTRParam[1]=svarg[20];
WTRParam[2]=svarg[21];
WTRParam[3]=svarg[22];
WTRParam[4]=0.0;

//cerr<<"Parameters: Calculate Stress:";
//for (int i=0; i<17; i++) cerr<<svarg[i]<<", ";
//cerr<<endl;


double Temp;

//cerr<<UI[0]<<' '<<UI[1]<<' '<<UI[2]<<' '<<UI[3]<<' '<<UI[4];
//this 2x3 lines are because the components in the added code are assumed in different sequence.
//probably this exchange is of no importance, but it is added for peace of mind

Temp=stress[4];
stress[4]=stress[5];
stress[5]=Temp;

Temp=D[4];
D[4]=D[5];
D[5]=Temp;

BBMPoint InitialPoint;
double StrainIncrement[7];

if (initial_p>0)
{
stress[0]=-initial_p+stress[0];
stress[1]=-initial_p+stress[1];
stress[2]=-initial_p+stress[2];
}


for (int i=0; i<6; i++)
{

    if (!isfinite(stress[i]))
    {
      cerr<<endl<<endl<<"Input stress ["<<i<<"] is not finite. Set to 0.001"<<endl<<endl<<endl;
      stress[i]=-0.001;
    }



    InitialPoint.stress[i]=-stress[i];
    StrainIncrement[i]=-D[i]*dt;
    if (!isfinite(PStar)) PStar=InitialPoint.GetMeanStress();
    if (!isfinite(suction)) suction=UI[13];
//    cerr<<"Stress:"<<-stress[0]<<' '<<-stress[1]<<' '<<-stress[2]<<' '<<-stress[3]<<' '<<-stress[4]<<' '<<-stress[5]<<endl;
//    cerr<<"Yield function:"<<MCModel.ComputeYieldFunctilon(&InitialPoint)<<endl;
//   cerr<<"Strain:"<<-StrainIncrement[0]<<' '<<-StrainIncrement[1]<<' '<<-StrainIncrement[2];
//    cerr<<' '<<-StrainIncrement[3]<<' '<<-StrainIncrement[4]<<' '<<-StrainIncrement[5]<<endl;
}


//Compute suction increment here, assuming undrained analysis and no water transport





    //cerr<<"Initial stress state"<<endl;
    //cerr<<"Stress:"<<-stress[0]<<' '<<-stress[1]<<' '<<-stress[2]<<' '<<-stress[3]<<' '<<-stress[4]<<' '<<-stress[5]<<endl;
    //cerr<<"Yield function:"<<MCModel.ComputeYieldFunctilon(&InitialPoint)<<endl;
   //cerr<<"Strain:"<<-StrainIncrement[0]<<' '<<-StrainIncrement[1]<<' '<<-StrainIncrement[2];
    //cerr<<' '<<-StrainIncrement[3]<<' '<<-StrainIncrement[4]<<' '<<-StrainIncrement[5]<<endl;

/*

  ps->getWithDefault("PMIN",UI[17],0.0);  //
  ps->getWithDefault("Accuracy",UI[18],0.0);  //
  ps->getWithDefault("CriticalSubstepSize",UI[19],0.0);  //
  ps->getWithDefault("DriftCorrection",UI[20],0.0);  //
  ps->getWithDefault("IntegrationAlgorithm",UI[21],0.0);  //

*/

//void BBMModel::SetModelParameters (double ShearModulusG, double Kappa_P, double Kappa_S, double PAtmospheric, double ReferenceStressPc,
//		double IncreaseCohesion_k, double ParameterDefLambS_r, double BetaDefLambdaS, double ParLambdaZero,
//		double ParNZero, double SlopeCSL_M, bool NonAssociatedFlowRule)
//cerr<<"G="<<G<<" K_p="<<KappaP<<" K_s="<<KappaS<<" p_atm="<<PAtmos<<" pc="<<pc<<" k="<<k<<" r="<<r<<" Beta="<<BetaBBM<<" Lambda="<<Lambda<<" N(0)="<<NZero<<" M="<<M<<" NonAssociated="<<NonAssociated<<endl;
    Model.SetDefaultIntegrationParameters();
    if (UI[17]>0) Model.PMIN=UI[17];
    if (UI[18]>0) Model.INTEGRATION_TOL=UI[18];
    if (UI[19]>0) Model.CriticalSubstepSize=UI[19];
    if (UI[20]>0) Model.DRIFT_CORRECTION=UI[20];
    if (UI[21]>0) Model.SOLUTION_ALGORITHM=UI[21];

    Model.SetModelParameters (G,KappaP, KappaS,PAtmos,pc,k,r,BetaBBM,Lambda,NZero,M,NonAssociated);
    InitialPoint.SetPStar(PStar);
    InitialPoint.SetSuction(suction);
    InitialPoint.SetSpecificVolume(Model.ComputeNu(InitialPoint.stress,InitialPoint.state,suction));


//cerr<<"Initial stress: ";
//for (int i=0; i<6; i++) cerr<<InitialPoint.stress[i]<<" ,";
//cerr<<endl;

//cerr<<"Pstar="<<InitialPoint.GetPStar()<<" suction="<<InitialPoint.GetSuction()<<" Nu="<<InitialPoint.GetSpecVol()<<endl;
//cerr<<"Strain increment:";
//for (int i=0; i<6; i++) cerr<<StrainIncrement[i]<<" ,";
//cerr<<endl;


    // BBMModel::SetModelParameters (double ShearModulusG, double Kappa_P, double Kappa_S, double PAtmospheric, double ReferenceStressPc,
	//	double IncreaseCohesion_k, double ParameterDefLambS_r, double BetaDefLambdaS, double ParLambdaZero,
	//	double ParNZero, double SlopeCSL_M, bool NonAssociatedFlowRule)



if (UseWaterRetention>0)
{
    double SpecVol=InitialPoint.GetSpecVol();
    double Sr=GetSr(UseWaterRetention,suction, WTRParam);   //get from current suction value
    double VolVoids=(1-Sr)*(SpecVol-1)/SpecVol;  //relative volume of voids
    double VolWater=Sr*(SpecVol-1)/SpecVol;
    double TotalVolumeInitial=VolVoids+VolWater;
    double dVol=StrainIncrement[0]+StrainIncrement[1]+StrainIncrement[2];
    double TotalVolume=TotalVolumeInitial-dVol;
    double SrNew=VolWater/TotalVolume;       //voids will reduce by dvol, so dSr=VolWater/TotalVolume
    if (SrNew>1.0) SrNew=1.0;
    double SuctionNew=GetSuction(UseWaterRetention,SrNew,WTRParam);                      //calculate from specified equation (van Genuchten or Gallipoli or...)
StrainIncrement[6]=SuctionNew-suction; //suction increment!!!
}
else StrainIncrement[6]=0.0; //suction increment!!!



int NoIter;

Model.Integrate (StrainIncrement,&InitialPoint,&NoIter);



//stress back for output
for (int i=0; i<6; i++)
{
        if (!isfinite(InitialPoint.stress[i]))
        {
          cerr<<endl<<endl<<"Output stress ["<<i<<"] is not finite. Set to 0.001"<<endl<<endl<<endl;
          InitialPoint.stress[i]=0.001;
        }
   // stress[i]=-(InitialPoint.stress[i]-initial_p);
    stress[i]=-InitialPoint.stress[i];
}

for (int i=0; i<3; i++)
{
    if (!isfinite(InitialPoint.state[i]))
    {
      cerr<<endl<<endl<<"Output state ["<<i<<"] is not finite"<<endl<<endl<<endl;
      getchar();
    }
}

//why those do not work???


/*
cerr<<"computed stress: ";
for (int i=0; i<6; i++) cerr<<stress[i]<<" ,";
cerr<<endl;
*/

//this 2x3 lines are because the components in the added code are assumed in different sequence.
//probably this exchange is of no importance, but it is added for the peace of mind

Temp=stress[4];
stress[4]=stress[5];
stress[5]=Temp;

Temp=D[4];
D[4]=D[5];
D[5]=Temp;


double MeanStress, LambdaS, Ns, SpecificVol, PZero;
MeanStress=InitialPoint.GetMeanStress();
if (InitialPoint.GetMeanStress()>InitialPoint.GetPStar()) InitialPoint.SetPStar(InitialPoint.GetMeanStress());
Ns=NZero-KappaS*log((suction+PAtmos)/PAtmos);
LambdaS=Lambda*((1-r)*exp(-BetaBBM*suction)+r);
PZero=(Lambda-KappaP)/(LambdaS-KappaP);
PZero=pc*pow( (InitialPoint.GetPStar() /pc),PZero);
SpecificVol=Ns-LambdaS*log(PZero/pc)+KappaP*log(PZero/MeanStress);




QuickK=InitialPoint.GetMeanStress()*SpecificVol/KappaP;
QuickNu=SpecificVol;

G=UI[0];
double Poisson=(3*QuickK-2*G)/(2*(3*QuickK+G));

if (Poisson>0.49)
{
  Poisson=0.49;
  G=3*QuickK*(1-2*Poisson)/(2*(1+Poisson));

}
else if (Poisson<0.01)
{
  Poisson=0.01;
  G=3*QuickK*(1-2*Poisson)/(2*(1+Poisson));

}



/*
for (int i=0; i<6; i++)
{
    if (fabs(InitialPoint.plastic_strain[i])>0.0)
        if (Factor<fabs(StrainIncrement[i]/InitialPoint.plastic_strain[i]))
                Factor=fabs(StrainIncrement[i]/InitialPoint.plastic_strain[i]);
}
*/

//copied from diamm

double Factor=1.0; //factor is  a quick fix, as otherwise the USM is too low and predicted stable time step is way too high

USM=Factor*(G+0.3*QuickK)/3.0;
//without Factor it
//seems that the USM is too low and the analysis fails. Not sure why
//as the elastic wave should be the quickest and it apparently work in diamm
//maybe I missed something important there


svarg[12]=InitialPoint.GetPStar();
svarg[13]=InitialPoint.GetSuction();
svarg[14]=0.0;
svarg[16]=QuickNu;
svarg[15]=QuickK;
svarg[17]=InitialPoint.GetMeanStress();
svarg[18]=InitialPoint.GetShearStress();
svarg[0]=G;

for (int i=0; i<19; i++)
{
    if (!isfinite(svarg[i]))
    {
      cerr<<endl<<endl<<"output svarg ["<<i<<"] is not finite"<<endl<<endl<<endl;
      cerr<<"Ns="<<Ns<<" Lambda="<<Lambda<<" LambdaS="<<LambdaS<<" Mean Stress="<<InitialPoint.GetMeanStress()<<" PStar="<<InitialPoint.GetPStar() <<" PZero="<<PZero<<" Nu="<<SpecificVol<<endl;
      getchar();
    }
}

//no neeed to update other svarg as those have not changed


}

double BBM::GetSr (double UseWaterRetention, double Suction, double * WTRParam)
{
    if (UseWaterRetention==1.0)
    {
        //VanGenuchten Model
        double m=WTRParam [0];
        double n=WTRParam [1];
        double alpha=WTRParam [2];
        double Sr=alpha*Suction;
        Sr=pow(Sr,n);
        Sr=1/(1+Sr);
        Sr=pow(Sr,m);
        if (Sr>1.0) Sr=1.0; if (Sr<0) Sr=0.0;
        return Sr;
    }
    else return 1.0;
}


double BBM::GetSuction (double UseWaterRetention, double Sr, double * WTRParam)
{
    if (UseWaterRetention==1.0)
    {
        //VanGenuchten Model
        double m=WTRParam [0];
        double n=WTRParam [1];
        double alpha=WTRParam [2];
        double Suction=pow(Sr,1/m);
        Suction= (1-Suction)/Suction;
        Suction=pow(Suction,1/n);
        Suction=Suction/alpha;
        return Suction;
    }
    else return 0;
}


void BBM::CheckModel(double UI[])
{
    //not implemented - some checks for input are in the constitutive_models.xml file

}
