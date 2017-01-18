using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;
using System.Text;

namespace ISAAR.MSolve.PreProcessor.Materials
{
    public class KavvadasClays : IFiniteElementMaterial3D
    {
        #region FieldsConstructor
        // Comments and explanations
        // Based on Kavvadas and Amorosi (2000) Belokas and Kalos investigations. We consider at first the edition of Kavvadas and Amorosi (2000) In the future it will be enriched with next editions. First Edition is on July 2016. Associative plasticity only in the current edition
        //Fields and properties
        public double Zeta { get; set; }
        public IMatrix2D<double> ConstitutiveMatrix { get; set; }
        public double[] Coordinates { get; set; }
        public int ID { get; set; }
        public bool Modified { get; set; }
        public double PoissonRatio { get; set; }
        public double[] Stresses { get; set; }
        public double YoungModulus { get; set; }
        public double ksi { get; set; }
        public double[] Qgrad = new double[6];
        public double[] Pgrad = new double[6];
        public double[] PAR;
        public double[] QH;
        public readonly double shearModulus;
        public readonly double[,] elasticConstitutiveMatrix;
        private double[] incrementalStrains = new double[6];
        private double plasticStrain;
        private double plasticStrainNew;
        private double[] stressesNew = new double[6];
        public void ClearState()
        {
            throw new NotImplementedException();
        }

        public void ClearStresses()
        {
            throw new NotImplementedException();
        }

        public object Clone()
        {
            throw new NotImplementedException();
        }

        public void ResetModified()
        {
            throw new NotImplementedException();
        }

        public void SaveState()
        {
            this.plasticStrain = this.plasticStrainNew;
            Array.Copy(this.stressesNew, this.Stresses, 6);
            //this.stresses = this.stressesNew.DeepClone();
        }

        public void UpdateMaterial(double[] strainsIncrement)
        {
            Array.Copy(strainsIncrement, this.incrementalStrains, 6);
            //this.incrementalStrains = strainsIncrement.DeepClone();
            //this.DecideLoad(PAR, QH, Stresses, elasticConstitutiveMatrix, ConstitutiveMatrix);
        }

        public KavvadasClays(double youngModulus, double poissonRatio, double alpha, double ksi)
        {
            this.YoungModulus = youngModulus;
            this.PoissonRatio = poissonRatio;
            //this.Alpha = alpha;
            this.ksi = ksi;
            this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
            double lamda = (youngModulus * poissonRatio) / ((1 + poissonRatio) * (1 - (2 * poissonRatio)));
            double mi = youngModulus / (2 * (1 + poissonRatio));
            double value1 = (2 * mi) + lamda;
            this.elasticConstitutiveMatrix = new double[6, 6];
            this.elasticConstitutiveMatrix[0, 0] = value1;
            this.elasticConstitutiveMatrix[0, 1] = lamda;
            this.elasticConstitutiveMatrix[0, 2] = lamda;
            this.elasticConstitutiveMatrix[1, 0] = lamda;
            this.elasticConstitutiveMatrix[1, 1] = value1;
            this.elasticConstitutiveMatrix[1, 2] = lamda;
            this.elasticConstitutiveMatrix[2, 0] = lamda;
            this.elasticConstitutiveMatrix[2, 1] = lamda;
            this.elasticConstitutiveMatrix[2, 2] = value1;
            this.elasticConstitutiveMatrix[3, 3] = mi;
            this.elasticConstitutiveMatrix[4, 4] = mi;
            this.elasticConstitutiveMatrix[5, 5] = mi;
        }

        public KavvadasClays(double youngModulus, double poissonRatio, double alpha, double ksi, double[] initialStresses) : this(youngModulus, poissonRatio, alpha, ksi)
        {
            var gamma = 20;
            var Kapa0 = 0.5;
            var Kmin = 0.1;
            var Kmax = 0.5;
            var Htot = 20;
            var Niso = 2.2;
            var sigmaref = 25.0;
            initialStresses[0] =-gamma*Zeta;
            initialStresses[1] =-gamma*Kapa0*Zeta;
            initialStresses[2] =-gamma*Kapa0*Zeta;
            initialStresses[3] =0;
            initialStresses[4] =0;
            initialStresses[5] =0;
            for (int i = 0; i < 6; i++)
                Stresses[i] = initialStresses[i];
            PAR[3] = (Kmin - Kmax) * Zeta / Htot + Kmax;
            var s0 = (Stresses[0] + Stresses[1] + Stresses[2]) / 3;
            QH[2] = Niso - PAR[3] * Math.Log(s0/sigmaref);
            PAR[0] = 0;
            PAR[1] = 0;
            PAR[2] = 0;
            PAR[4] = 0;
            PAR[5] = 0;
            PAR[6] = 0;
            PAR[7] = 0;
            PAR[8] = 0;
            PAR[9] = 0;
            PAR[10] = 0;
            PAR[11] = 0;
            PAR[12] = 0;
            PAR[13] = 0;
            PAR[14] = 0;
            PAR[15] = 0;
            PAR[16] = 0;
            PAR[17] = 0;
            PAR[18] = 0;
            PAR[19] = 0;
            QH[0] = 0;
            QH[1] = 0;
            QH[3] = 0;
            QH[4] = 0;
            QH[5] = 0;
            QH[6] = 0;
            QH[7] = 0;
            QH[8] = 0;
            QH[9] = 0;
            QH[10] = 0;
            QH[11] = 0;
            QH[12] = 0;
            QH[13] = 0;
            QH[14] = 0;
            QH[15] = 0;
            QH[16] = 0;
            QH[17] = 0;
            QH[18] = 0;
            QH[19] = 0;
            QH[20] = 0;
            //Dont forget that PAR[0] and QH[0] here are not used (see original code)
        }

        #endregion
            #region OriginalCode
            //        C CYCLIC.FOR
            //C
            //C   Revision 30/11/2012: by Al.Kalos: MCC W/ PYE + ISE + DESTRUCTURING implementation
            //C Revision (resumed) 30.01.2013 (Al.Kalos)
            //C
            //C********************************************************************************************
            //C   Subroutine library for the new cyclic model
            //C Works with transformed stresses - strains
            //C   Compressive stresses - strains are positive
            //C********************************************************************************************
            //C********************************************************************************************
            //C   Main driver routines:
            //C CYCLIC : Finite strain driver.Given DE updates S, QH and computes C.
            //C SUBROUTINE CYCLIC(DLIMIT, NT, PAR, DE, QH, S, C)
            //C CSTRAIN: Infinitesmal strain driver.Given DE updates S, QH.
            //C called by CYCLIC.
            //C                 SUBROUTINE CSTRAIN (NT, PAR, DE, QH, S)
            //C CSTIF  : Infinitesmal strain driver.Computes appropriate(elastic or
            //C elastoplastic) stiffness C.Called by CYCLIC.
            //C SUBROUTINE CSTIF (ISTIF, NT, QH, S, PAR, C)
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE CYCLIC(DLIMIT, NT, PAR, DE, QH, S, C)
            //C Finite-strain driver for incremental models.
            //C   Given the current state (S, QH) and an applied finite
            //C strain increment(DE) (all in the transformed space, with soil mechanics
            //C   sign conventions), the routine:
            //C Updates the state(S, QH) and computes the average continuum Jacobian C.
            //C The routine splits a finite strain increment DE in a number NDIV of
            //C infinitesmal strain increments DDE(DDE<DLIMIT) and repeatedly applies
            //C   all of them.
            //C
            //C   Input:  DLIMIT,NT,PAR,DE,QH,S
            //C   Output: QH,S,C
            //C
            //C DLIMIT: Max value of the strain increment(max.component) to be
            //C considered infinitesmal.If any of the applied DE(I) exceeds
            //C           DLIMIT, the strain increment is sub-divided.
            //C NT :    Size of the stress space(2=triaxial test, 3=plane strain test,
            //C           4=DSS test, 4=plane strain and axisymmetric problems,
            //C           6=general 3-D problems
            //C
            //C   PAR(dim= NT + 11):    stores material constants
            //C   PAR(1)=            In poro-elasticity: (2G/K)
            //C PAR(2)=            Lamdastar - slope of the ICL(Cam-Clay compression)
            //C PAR(3)=            In poro-elasticity: kappastar(Cam-Clay)
            //C PAR(4)...(NT+2)=   c1, c2, ... c(NT-1)=ratio of ellipse axes
            //C PAR(NT+3)=B0 B0>=Bres>=1   (ratio of initial structure)  B0=a/a*
            //C   PAR(NT+4)=Bres Bres>=1       (ratio of residual structure) Bres=a_res/a*
            //C   PAR(NT+5)=n_vp(volumetric destructuring variable)
            //C PAR(NT+6)=n_qp(deviatoric destructuring variable)
            //C PAR(NT+7)=z_vp(volumetric overshooting variable)
            //C PAR(NT+8)=z_qp(deviatoric overshooting variable)
            //C PAR(NT+9)=ksi(ratio of SSE/PYE ellipsoidal shapes)
            //C PAR(NT+10)=L1(constant employed in the interpolation for the PYE plastic modulus)
            //C PAR(NT+11)=gamma(exponent employed in the interpolation for the PYE plastic modulus)
            //C PAR(NT+12)=OCR OCR
            //C PAR(NT+13)=Niso
            //C
            //C DE(dim= NT):        Finite strain increment in transformed space
            //C
            //C   QH(dim= 2 * NT + 8):    stores all state variables(they vary during deformation)
            //C QH(1)=             alpha - Center of the Yield Surface(SSE)
            //C QH(2)=             specific volume(for poro-elastic volumetric hardening)
            //C QH(3)=             NSURF(0=inside YS, -1=on PYE, +1=on SSE)
            //C QH(4)=             plastic volumetric strain
            //C   QH(5)=             plastic deviatoric strain
            //C   QH(6)=             alphastar - Center of the Intrinsic Strength Envelope(ISE==SSE)
            //C QH(7)...QH(NT+6)=  initial state of stress on PYE(S) once I enter the plastic region
            //C QH(NT+7)...QH(2*NT+6)=center of yield surface PYE(L)
            //C QH(NT+7)...QH(2*NT+6)=center of yield surface PYE(L)
            //C QH(2*NT+7)         H'' at the conjugate point
            //C QH(2*NT+8)         H_partnew -> H=H''+H_partnew
            //C
            //C S(dim= NT) =         current stress in transformed space: sigma, S1, S2, S(NT-1)
            //C C(dim= NT * NT) =      continuum Jacobian
            //C----------------------------------------------------------------
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION QH(1),PAR(1),S(1),DE(1)
            //      DIMENSION DDE(6),C(36),CC(36),SS(6)
            //C
            //C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            //C********************************************************************************************
            //C*********************************** USER INTERFACE*************************************
            //C********************************************************************************************
            //C   Note: Size of array QQH must be sufficient to store array QH,
            //C Currently = 2 * NT + 8 where max of NT = 6
            //C
            //      DIMENSION QQH(2*NT+8)
            //      NQH=2*NT+8
            //C
            //C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            //C
            //      NN = NT * NT
            //C Zero C matrix
            //      DO 5 I=1,NN
            //    5 C(I)=0.D0
            //C   Compute NDIV and inf.strain subincrement DDE
            //      DMAX= 0.D0
            //      DO 10 I= 1, NT
            //      DUM = DABS(DE(I))
            //      IF (DUM.GT.DMAX) DMAX=DUM
            //   10 CONTINUE
            //      NDIV = 1
            //      IF(DMAX.GT.DLIMIT) NDIV=IDINT(DMAX/DLIMIT+1.D0)
            //      DDIV=DBLE(NDIV)
            //      DO 20 I=1,NT
            //   20 DDE(I)=DE(I)/DDIV
            //C   Repeatedly apply DDE, NDIV times
            //C  --------------------------------------------------------------------------------------------
            //      DO 100 INC=1,NDIV
            //C   Save input state:
            //      DO 35 I=1,NT
            //   35 SS(I)=S(I)
            //      DO 36 I=1,NQH
            //   36 QQH(I)=QH(I)
            //C Apply DDE and compute new state(QH, S):
            //      CALL CSTRAIN(NT, PAR, DDE, QH, S)
            //C Compute continuum Jacobian CC at initial state(QQH, SS).
            //C Check if elastic or elastoplastic stiffness should be computed.
            //    IF (QH(3).EQ.0.D0) THEN
            //        ISTIF = 0
            //      ELSE
            //        IF(QQH(3).EQ.0.D0) THEN
            //          ISTIF = 0
            //        ELSE
            //           ISTIF = 1
            //        END IF
            //      END IF
            //      CALL CSTIF(ISTIF, NT, QQH, SS, PAR, CC)
            //C Add effect of Jacobian CC to running total C:
            //      DO 50 I=1,NN
            //   50 C(I)=C(I)+CC(I)
            //  100 CONTINUE
            //C  --------------------------------------------------------------------------------------------
            //C Compute mean stiffness
            //      DO 120 I=1,NN
            //  120 C(I)=C(I)/DDIV
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE CSTRAIN(NT, PAR, DE, QH, S)
            //C Given the current state(S, QH) and an applied infinitesmal
            //C strain increment(DE) (all in the transformed space), the routine
            //C updates the state(S, QH).
            //C Input:  NT,PAR,DE,S,QH
            //C   Output:  S,QH
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION QH(1),PAR(1),S(1),DE(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Check if current state is inside the yield surface(NSURF= 0)
            //      IF(QH(3).EQ.0.D0) THEN
            //C   Current state is inside the yield surface:
            //         CALL UNLOAD(NT, PAR, DE, QH, S)
            //      ELSE
            //C   Current state is on the yield surface(NSURF= -1 or +1).
            //C Check if the strain increment causes loading or unloading
            //C   by computing XNUM=Q.C.DE
            //       CALL NUMER(NT, PAR, DE, QH, S, XNUM)
            //         IF(XNUM.LT.0.D0) THEN
            //C   Unload
            //           CALL UNLOAD(NT, PAR, DE, QH, S)
            //         ELSEIF(XNUM.EQ.0.D0) THEN
            //C   In case of neutral loading call LOADNEUTRAL
            //           CALL LOADNEUTRAL(NT, PAR, DE, QH, S)
            //         ELSE
            //            CALL LOAD(NT, PAR, DE, QH, S)
            //         END IF
            //      END IF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE LOAD(NT, PAR, DE, QH, S)
            //C
            //C   Given an initial state(S, QH) on the yield surface(and either on or
            //C inside the bounding surface) and an applied
            //C   infinitesmal strain increment(DE) -which is known to cause plastic
            //C loading-, the routine updates the state(S, QH).
            //C Note: on input: QH(NT+3)=NSURF= +1 or -1.
            //C Input: NT,PAR,DE,QH,S
            //C   Output: QH,S
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION QH(1),PAR(1),S(1),DE(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //      IF(QH(3).EQ.-1.D0) THEN
            //C   Stress point is inside the SSE(but on PYE):
            //         CALL YLOAD(NT, PAR, DE, QH, S)
            //      ELSE
            //C   Stress point is simultaneously on the SSE and PYE:
            //         CALL BLOAD(NT, PAR, DE, QH, S)
            //      END IF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE BLOAD(NT, PAR, DE, QH, S)
            //C The current stress state S lies simultaneously
            //C on the PYE and the SSE
            //C   The strain increment DE causes plastic loading.
            //C The routine updates state(QH, S).
            //C Input: NT,PAR,DE,QH,S
            //C   Output: QH,S
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION QH(1),PAR(1),S(1),DE(1)
            //      DIMENSION Q(6), P(6), DD(6), DS(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute XLDOT=Lambda-dot
            //    CALL COMPLDOT(NT, PAR, QH, S, DE, XLDOT)
            //C Compute the Pf of the PYE(S on the SSE)
            //      CALL COMPP(NT, PAR, QH, S, P)
            //C Compute the DSQRT(2/3*P':P')
            //      CALL COMPDOTPDOT(NT, P, DPP)
            //C Compute DEvp and DEqp
            //    DEVP = 0.D0
            //      DEQP=0.D0
            //      DEVP = DABS(XLDOT * P(1))
            //      DEQP=DABS(XLDOT* DPP)
            //C Compute DEvp_true=XLDOT* P(used for hardening or softening -> NOT DEST)
            //      DEVPTRUE=0.D0
            //      DEVPTRUE = DABS(XLDOT) * P(1)
            //C Compute Evp and Eqp
            //    COMPEVP = QH(4) + DEVP
            //      COMPEQP=QH(5)+DEQP
            //C   DEXPART=EXP(-(n_vp* E_vp+n_qp* E_qp))
            //      DEXPART=DEXP(-DABS(PAR(NT+5)*COMPEVP+PAR(NT+6)*COMPEQP))
            //C DPART1 = B0 - Bres + z_vp * E_vp + z_qp * E_qp
            //      DPART1=
            //     1 DABS(PAR(NT+3)-PAR(NT+4)+PAR(NT+7)*COMPEVP+PAR(NT+8)*COMPEQP)
            //C  -------------------------------------------------------------------------------------------- 
            //C Update Plastic Volumetric Strain
            //    QH(4)=COMPEVP
            //C   Update Plastic Deviatoric Strain
            //      QH(5)=COMPEQP
            //C   Compute ALPHA of the ISE bubble -> use void ratio of previous inc. -- not the updated QH(2)
            //      QH(6)=QH(6)*(1.D0+QH(2)*DEVPTRUE/(PAR(2)-PAR(3))) 
            //C Update center of bounding surface
            //      QH(1)=QH(6)*(DPART1* DEXPART+PAR(NT+4))
            //C Compute elastic strain increment DD
            //      DO 27 I=1,NT
            //   27 DD(I)=DE(I)-DABS(XLDOT)*P(I)
            //C Compute stress increment DS for DD
            //    CALL COMPDS(NT, PAR, QH, S, DD, DS)
            //C Update stress & center of PYE
            //    DO 30 I=1,NT
            //C   Update stress
            //   30 S(I)=S(I)+DS(I)
            //C Adjust stress(S) on the bounding surface
            //      ALPHA=QH(1)
            //C Update specific volume
            //      QH(2)=QH(2)*(1.D0-DE(1))
            //C  --------------------------------------------------------------------------------------------
            //C I will run the ADJUST
            //C  --------------------------------------------------------------------------------------------
            //C Adjust stress(S) on the SSE
            //   CALL ADJBSE(NT, PAR, ALPHA, S)
            //C Compute center of yield surface
            //      CALL COMPSL(NT, PAR, S, QH)
            //C NSURF was 1 and preserves its value
            //      QH(3)=1.D0
            //C  --------------------------------------------------------------------------------------------
            //C The main text of the subroutine has finished -> I will run the check for future debugging purposes
            //C   Compute the Fvalue of SSE(=FLAG)
            //      CALL FBSE(NT, S, PAR, QH(1),FLAG)
            //C Check whether FLAG is a real number or NaN
            //    IF(DABS(FLAG).GE.0.D0) THEN
            //   GOTO 28
            //      ELSE
            //      WRITE(6,*) 'BLOAD SUBROUTINE FAILED LINE 604'
            //        STOP 2
            //      ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //   28 RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE YLOAD(NT, PAR, DE, QH, S)
            //C The current stress state S lies on the PYE and
            //C   the center of the PYE(SL) inside the bounding surface.
            //C   Strain increment DE causes plastic loading.
            //C The routine updates state (S, QH).
            //C Input NT,PAR,DE,QH,S
            //C   Output: QH,S
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION QH(1), PAR(1), S(1), DE(1)
            //      DIMENSION P(6), B(6), DS(6), DE2(6), SN(6), DE1(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute XLDOT=Lambda-dot.Modulus H is computed by using an
            //C interpolation from the conjugate stress point on the SSE.
            //    CALL COMPLDOT (NT, PAR, QH, S, DE, XLDOT)
            //C Compute array Pf
            //      CALL COMPP(NT, PAR, QH, S, P)
            //C Compute the deviatoric part of Q
            //    CALL COMPDOTPDOT(NT, P, DPP)
            //C Compute elastic strain increment DE2
            //      DO 27 I=1,NT
            //   27 DE2(I)=DE(I)-DABS(XLDOT)* P(I)
            //C Compute stress increment DS for DE2
            //    CALL COMPDS(NT, PAR, QH, S, DE2, DS)
            //C  --------------------------------------------------------------------------------------------
            //C Compute B(direction of translation of yield surface)
            //C  --------------------------------------------------------------------------------------------
            //C Compute volumetric component of B
            //      B(1)=(S(1)-QH(NT+7))/PAR(NT+9)-(S(1)-QH(1))
            //      DO 40 I=2,NT
            //C   Compute deviatoric component of B
            //   40 B(I)=(S(I)-QH(NT+I+6))/PAR(NT+9)-S(I)
            //C  --------------------------------------------------------------------------------------------
            //C Compute DEvp and DEqp
            //    DEVP = DABS(XLDOT * P(1))
            //      DEQP=DABS(XLDOT* DPP)
            //C Compute DEvp_true=XLDOT* P(used for hardening or softening -> NOT DEST)
            //      DEVPTRUE=DABS(XLDOT)* P(1)
            //C Compute Evp and Eqp
            //    COMPEVP = QH(4) + DEVP
            //      COMPEQP=QH(5)+DEQP
            //C   DEXPART=EXP(-(n_vp* E_vp+n_qp* E_qp))
            //      DEXPART=DEXP(-DABS(PAR(NT+5)* COMPEVP+PAR(NT+6)* COMPEQP))
            //C DPART1 = B0 - Bres + z_vp * E_vp + z_qp * E_qp
            //      DPART1=PAR(NT+3)-PAR(NT+4)+PAR(NT+7)* COMPEVP+PAR(NT+8)* COMPEQP
            //      DPART1=DABS(DPART1)
            //C  -------------------------------------------------------------------------------------------- 
            //C Compute ALPHA of the ISE bubble -> use void ratio of previous inc. -- not the updated QH(2)
            //      QH6=QH(6)*(1.D0+QH(2)* DEVPTRUE/(PAR(2)-PAR(3)))
            //C Update center of bounding surface
            //      ALPHA=QH6*(DPART1* DEXPART+PAR(NT+4))
            //C Compute the DQH(1)=DA
            //    DA = ALPHA - QH(1)
            //C
            //C  --------------------------------------------------------------------------------------------
            //C Compute the ratio a_dot/a
            //C  --------------------------------------------------------------------------------------------
            //      RATIO=DA/QH(1)
            //C  --------------------------------------------------------------------------------------------
            //C Compute YPLDOT=1.D0+a_dot/a=a_NEW/a
            //    YPLDOT = 1.D0 + RATIO
            //      DO 10 I=1,NT
            //C   Update stress S
            //   10 SN(I)=S(I)+DS(I)
            //C Check the Fvalue of the SSE(=F_CHECKPOINT)
            //      CALL FBSE(NT, SN, PAR, ALPHA, F_CHECKPOINT)
            //C Compute YLDOT - DM:
            //      CALL COMPYLDOT(NT, PAR, S, QH, DS, DE, YLDOT, DA)
            //      YPLDOT=1.D0+RATIO
            //C  --------------------------------------------------------------------------------------------
            //C Check whether l_dot>=1.D0+a_dot/a -> If so it means that I have exceeded the
            //C   the portion of B, in which case I have crossed the SSE
            //C  --------------------------------------------------------------------------------------------
            //C
            //       IF(F_CHECKPOINT.LT.0.D0) THEN
            //C
            //C  --------------------------------------------------------------------------------------------
            //C The PYE remains inside the SSE
            //C
            //C Update the size of SSE
            //       QH(1)=ALPHA
            //C
            //       DO 50 I=1,NT
            //C   Update center of PYE
            //       QH(NT+I+6)=YPLDOT* QH(NT+I+6)+YLDOT* B(I)
            //C Update stress S
            //   50  S(I)=S(I)+DS(I)
            //C Compute Evp and Eqp
            //     QH(4)=QH(4)+DEVP
            //     QH(5)=QH(5)+DEQP
            //C   Compute ALPHA of the ISE bubble
            //       QH(6)=QH6
            //C   Update specific volume
            //       QH(2)=QH(2)*(1-DE(1))
            //C  --------------------------------------------------------------------------------------------
            //C Adjust stress(S) on the bounding surface
            //       CALL ADJPYE(NT, PAR, QH, S)
            //C Update NSURF
            //QH(3)=-1.D0
            //C  --------------------------------------------------------------------------------------------
            //      ELSE
            //C  --------------------------------------------------------------------------------------------
            //C The PYE tends to cross the SSE.
            //C   Limit the strain increment.
            //C Compute fraction RR(<=1) of RATIO which just reaches the SSE.
            //C
            //     CALL INTBSE(NT, S, PAR, QH, DS, RR)
            //C
            //       DO 61 I=1,NT
            //       DE2(I)=(1.D0-RR)* DE(I)
            //   61  DE1(I)=RR* DE(I)
            //C Compute XLDOT=Lambda-dot.Modulus H is computed by using an
            //C interpolation from the conjugate stress point on the SSE.
            //    CALL COMPLDOT (NT, PAR, QH, S, DE1, XLDOT)
            //C Compute array Pf
            //      CALL COMPP(NT, PAR, QH, S, P)
            //C Compute the deviatoric part of Q
            //    CALL COMPDOTPDOT(NT, P, DPP)
            //C Compute DEvp and DEqp
            //    DEVP = DABS(XLDOT * P(1))
            //      DEQP=DABS(XLDOT* DPP)
            //C Compute DEvp_true=XLDOT* P(used for hardening or softening -> NOT DEST)
            //      DEVPTRUE=DABS(XLDOT)* P(1)
            //C Compute Evp and Eqp
            //    COMPEVP = QH(4) + DEVP
            //      COMPEQP=QH(5)+DEQP
            //C   DEXPART=EXP(-(n_vp* E_vp+n_qp* E_qp))
            //      DEXPART=DEXP(-DABS(PAR(NT+5)* COMPEVP+PAR(NT+6)* COMPEQP))
            //C DPART1 = B0 - Bres + z_vp * E_vp + z_qp * E_qp
            //      DPART1=PAR(NT+3)-PAR(NT+4)+PAR(NT+7)* COMPEVP+PAR(NT+8)* COMPEQP
            //      DPART1=DABS(DPART1)
            //C  -------------------------------------------------------------------------------------------- 
            //C Compute ALPHA of the ISE bubble -> use void ratio of previous inc. -- not the updated QH(2)
            //      QH(6)=QH(6)*(1.D0+QH(2)* DEVPTRUE/(PAR(2)-PAR(3)))
            //C Update center of bounding surface
            //      QH(1)=QH(6)*(DPART1* DEXPART+PAR(NT+4))
            //C Compute Evp and Eqp
            //     QH(4)=QH(4)+DEVP
            //     QH(5)=QH(5)+DEQP
            //C   Update specific volume
            //     QH(2)=QH(2)*(1-DE1(1))
            //C Update stress S
            //       DO 60 I=1,NT
            //   60  S(I)=S(I)+RR* DS(I)
            //C  --------------------------------------------------------------------------------------------
            //C Adjust stress(S) on the bounding surface
            //       CALL ADJBSE(NT, PAR, QH(1),S)
            //C Place center of yield surface so that yield surface touches SSE at S
            //       CALL COMPSL(NT, PAR, S, QH)
            //C Update NSURF:
            //       QH(3)=1.D0
            //C   If there is a portion of the strain increment left, apply it.
            //       IF(RR.LT.1.D0) THEN
            //C  --------------------------------------------------------------------------------------------
            //C Aply remaining increment
            //C  --------------------------------------------------------------------------------------------
            //C Impose the remaining strain increment using BLOAD
            //         CALL BLOAD(NT, PAR, DE2, QH, S)
            //C  --------------------------------------------------------------------------------------------
            //       ENDIF
            //C  --------------------------------------------------------------------------------------------
            //      ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C The main text of the subroutine has finished -> I will run the check for future debugging purposes
            //C   Compute the Fvalue of PYE(=FLAG)
            //      CALL FPYE(NT, S, PAR, QH, FLAG)
            //C Check whether FLAG is a real number or NaN
            //    IF(DABS(FLAG).GE.0.D0) THEN
            //   GOTO 28
            //      ELSE
            //      WRITE(6,*) 'YLOAD SUBROUTINE FAILED LINE 785'
            //        STOP 2
            //      ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //   28 RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE UNLOAD(NT, PAR, DE, QH, S)
            //C Routine performs unloading(i.e.movement of the stress point towards
            //C the interior of the yield surface), starting from a state on the yield
            //C   surface or inside the yield surface.
            //C   If the new computed state is inside the yield surface, deformation is
            //C purely elastic.
            //C If the new computed state has crossed the yield surface, deformation is
            //C initially elastic and then plastic.In this case, the routine computes
            //C   the point of intersection (SI) with the yield surface,
            //C   and performs plastic loading for the remainder of the strain increment.
            //C   Returns the updated state: (QH, S).
            //C Input:  NT,PAR,DE,QH,S
            //C   Output: QH,S
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1), QH(1), DE(1), S(1)
            //      DIMENSION DE2(6), DS(6), SN(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute stress increment(DS), and increment of the specific volume(DV),
            //C assuming that deformation(DE) is purely elastic.
            //   CALL COMPDS (NT, PAR, QH, S, DE, DS)
            //      DV=-DE(1)* QH(2)
            //C Create an updated(test) stress SN but do not store yet -> CHECK first
            //      DO 10 I=1,NT
            //   10 SN(I)=S(I)+DS(I)
            //C Check if new stress state(SN) is inside the yield surface
            //      ALPHA=QH(1)
            //C Compute the Fvalue of PYE(=F_CHECKPOINT2)
            //      CALL FPYE(NT, SN, PAR, QH, F_CHECKPOINT2)
            //C  --------------------------------------------------------------------------------------------
            //      IF(F_CHECKPOINT2.LT.0.D0) THEN
            //C  --------------------------------------------------------------------------------------------
            //C Final state is inside the yield surface.
            //C   Update stress and hardening parameters and return

            //      DO 20 I= 1, NT
            //   20   S(I)=SN(I)
            //C Update specific volume
            //        QH(2)=QH(2)+DV
            //C   Update NSURF
            //        QH(3)=0.D0
            //C  --------------------------------------------------------------------------------------------
            //      ELSE
            //C  --------------------------------------------------------------------------------------------
            //C Final state is outside or on the yield surface.
            //C   Compute ratio XL which defines the point SI on the PYE
            //C   in the direction of the stress increment DS:
            //        CALL INTPYE (NT, S, PAR, QH, DS, XL3)
            //C Update stress state and compute remaining strain increment(DE2):
            //        DO 30 I=1,NT
            //C   Update Stresses
            //        S(I)=S(I)+XL3* DS(I)
            //C Compute remaining increment
            //   30   DE2(I)=(1.D0-XL3)* DE(I)
            //C Update specific volume
            //        QH(2)=QH(2)+XL3* DV
            //C Adjust State of Stress S on the PYE
            //      CALL ADJPYE(NT, PAR, QH, S)
            //C Call the FVALUE to check whether the center of PYE lays on the BSE
            //      CALL FBSE(NT, S, PAR, ALPHA, F_CHECKPOINT3)
            //C  --------------------------------------------------------------------------------------------
            //C Check whether the center of PYE lays on the BSE
            //C  --------------------------------------------------------------------------------------------
            //        IF(F_CHECKPOINT3.LT.0.D0) THEN
            //         QH(3)=-1.D0
            //       ELSE
            //C Adjust stress(S) on the bounding surface
            //          CALL ADJBSE(NT, PAR, ALPHA, S)
            //C Place center of yield surface so that yield surface touches SSE at S
            //          CALL COMPSL(NT, PAR, S, QH)
            //C Update NSURF
            //        QH(3)=1.D0
            //      ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C Store the Initial Stress S -> entered plastic zone
            //      DO 50 I=1,NT
            //   50   QH(I+6)=S(I)
            //C  --------------------------------------------------------------------------------------------
            //C Apply remaining strain increment DE2: compute final state.
            //      CALL LOAD(NT, PAR, DE2, QH, S)
            //C  --------------------------------------------------------------------------------------------
            //      END IF
            //C  --------------------------------------------------------------------------------------------
            //C The main text of the subroutine has finished -> I will run the check for future debugging purposes
            //C   Compute the Fvalue of PYE(=FLAG)
            //      CALL FPYE(NT, S, PAR, QH, FLAG)
            //C Check whether FLAG is a real number or NaN
            //    IF(DABS(FLAG).GE.0.D0) THEN
            //   GOTO 28
            //      ELSE
            //      WRITE(6,*) 'UNLOAD SUBROUTINE FAILED LINE 886'
            //        STOP 2
            //      ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //   28 RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE LOADNEUTRAL(NT, PAR, DE, QH, S)
            //C Routine performs neutral loading(i.e.movement of the stress point vertical
            //C to the gradient of the yield surface), starting from a state on the yield
            //C   surface.
            //C The new state(same as the old) is on the yield surface, but the deformation is
            //C purely elastic.
            //C Returns the updated state: (QH,S).
            //C Input:  NT,PAR,DE,QH,S
            //C   Output: QH,S
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1), QH(1), DE(1), S(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C The increment of the specific volume(DV)
            //      DV=-DE(1)* QH(2)
            //C Update specific volume
            //        QH(2)=QH(2)+DV
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE NUMER(NT, PAR, DE, QH, S, XNUM)
            //C Computes the numerator of Lambda-dot= XNUM=Q* C* DE
            //C Input:  NT,PAR,DE,QH,S
            //C   Output: XNUM
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1),QH(1),DE(1),S(1)
            //      DIMENSION Q(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute vector Q wrt the yield surface
            //      CALL COMPQPYE(NT, PAR, QH, S, Q)
            //C
            //      CALL QUAD(NT, PAR, QH, S, Q, DE, XNUM)
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE FBSE(NT, S, PAR, ALPHA, F)
            //C Computes the value(FBSE) of the bounding function at the
            //C state(S). Surface defined by PAR, ALPHA
            //C Input: NT,S,PAR,ALPHA
            //C   Output: FVALUE
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION S(1), PAR(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //      F=0.D0
            //C  --------------------------------------------------------------------------------------------
            //      DO 10 I=2,NT
            //   10 F=F+DABS(S(I)/PAR(I+2))* DABS(S(I)/PAR(I+2))
            //      F=F+DABS((ALPHA-S(1))*(ALPHA-S(1)))-DABS(ALPHA* ALPHA)
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE FPYE(NT, S, PAR, QH, F)
            //C Computes the value(F) of the yield/bounding function at the
            //C state(S). Surface defined by PAR, ALPHA
            //C Input: NT,S,PAR,QH
            //C   Output: F
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION S(1), PAR(1), QH(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //      F=DABS((S(1)-QH(NT+7))*(S(1)-QH(NT+7)))
            //      DHALFTERM=0.D0
            //      DO 10 I=2,NT
            //      DHALFTERM = (S(I) - QH(NT + I + 6)) / PAR(I + 2)
            //   10 F=F+DABS(DHALFTERM* DHALFTERM)
            //      F=F-DABS((PAR(NT+9)* QH(1))*(PAR(NT+9)* QH(1)))
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE ADJBSE(NT, PAR, ALPHA, S)
            //C Given a stress state(S) inside, which is close to(but maybe not
            //C   exactly on the bounding surface, due to the finiteness of the strain
            //C   increment) the routine updates the stress state(S) and places it on
            //C the surface, by projecting it with pole the center of the surface
            //C
            //C   Input: NT,PAR,ALPHA
            //C   Output: S
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION S(1), PAR(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C FVALUE = F(SIGMA, ALPHA)
            //      CALL FBSE(NT, S, PAR, ALPHA, FVALUE)
            //C XLAMDA = DSQRT(ALPHA * ALPHA / (FVALUE + ALPHA * ALPHA))
            //      ALPHA2=DABS(ALPHA* ALPHA)
            //      XLAMDA=DSQRT(DABS(ALPHA2/(FVALUE+ALPHA2)))
            //C Check whether XLAMDA is real or NaN
            //    IF(DABS(XLAMDA).GE.0.D0) THEN
            //C    Adjust stress state on bounding surface
            //      S(1)=ALPHA+XLAMDA*(S(1)-ALPHA)
            //      DO 30 I=2,NT
            //   30 S(I)=XLAMDA* S(I)
            //      ELSE
            //        WRITE(6,*) 'ADJBSE LINE 1008 FUCKED'
            //        STOP 11
            //      ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE ADJPYE(NT, PAR, QH, S)
            //C Given the current stress state(S), which is close to(but maybe
            //C not exactly on the surface, due to the finiteness of the strain
            //C   increment) the routine updates the stress state(S) and places it
            //C   on the surface, by projecting it with pole the center of the surface.
            //C   Input: NT, PAR, QH, S
            //C   Output: S
            //C

            //  IMPLICIT REAL*8 (A-H, O-Z)
            //      DIMENSION S(1), PAR(1), QH(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //      CALL FPYE(NT, S, PAR, QH, FVALUE)
            //C  --------------------------------------------------------------------------------------------
            //      DUM=DABS((PAR(NT+9)* QH(1))*(PAR(NT+9)* QH(1)))
            //C  --------------------------------------------------------------------------------------------
            //      XD=DUM+FVALUE
            //      IF(FVALUE.NE.0.D0) THEN
            //       XL2 = DSQRT(DABS(DUM / XD))
            //C Check whether XL2 is a real number or NaN
            //      IF(DABS(XL2).GE.0.D0) THEN
            //       DO 11 I=1,NT
            //   11     S(I)=XL2* S(I)+(1.D0-XL2)* QH(NT+I+6)
            //        ELSE
            //          WRITE(6,*) 'ADJPYE LINE 1045 FUCKED'
            //          STOP 12
            //        ENDIF
            //      ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE INTBSE(NT, S, PAR, QH, DS, XL)
            //C Given a stress state(S) inside the BSE(on the PYE) and a stress increment
            //C(DS) such that(S+DS) is outside the BSE, the routine computes
            //C   the factor XL(0..1) where the stress S=S+XL*DS, lays on the BSE (and on the PYE).
            //C   Also, given a stress state (S) on the BSE (i.e.CC=0),
            //C   and a stress increment (DS) such that (S+DS) is outside the BSE,
            //C   but gone through unloading (i.e.crossing the inside of the PYE),
            //C   the routine computes the factor XL (0..1) where the stress S=S+XL*DS,
            //C   is on the PYE surface.
            //C
            //C   Input: NT, S, PAR, QH, DS
            //C   Output: XL
            //C

            //  IMPLICIT REAL*8 (A-H, O-Z)

            //  DIMENSION S(1), DS(1), PAR(1), QH(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C   Initialize AA and BB

            //  AA=DABS(DS(1)*DS(1))

            //  BB=DS(1)*(S(1)-QH(1))
            //C   Finalize AA and BB

            //  DO 10 I=2, NT

            //  BB=BB+((DS(I)*S(I))/(PAR(I+2)*PAR(I+2)))
            //   10 AA=AA+DABS(DS(I)/PAR(I+2))*DABS(DS(I)/PAR(I+2))
            //C  --------------------------------------------------------------------------------------------
            //C   Check AA

            //  IF (AA.EQ.0.D0) THEN

            //      WRITE (6,101)
            //          WRITE (11,101)
            //  101     FORMAT (1X,'SUBROUTINE INTERSECT INTBSE: VECTOR DS=0')
            //          XL=1.D0

            //      GOTO 28

            //  END IF
            //C  --------------------------------------------------------------------------------------------
            //C   Compute the Fvalue of SSE (=CC)

            //  CALL FBSE (NT, S, PAR, QH(1), CC)
            //C   Compute SUBSQUARE

            //  SUBSQUARE=DABS(BB*BB)-AA*CC
            //C   Compute XL

            //  XL=(-BB+DSQRT(DABS(SUBSQUARE)))/AA
            //C   Check if XL is a real number or NaN

            //  IF (DABS(XL).GE.0.D0) THEN
            //C   If F(s)=CC=0 we want the non-zero root:

            //    IF (CC.EQ.0.D0) XL=1.D0

            //    IF (XL.GE.1.D0) XL=1.D0

            //    IF (XL.LT.0.D0) XL=0.D0

            //  ELSE

            //    WRITE (6,*) 'INTBSE LINE 1102 FUCKED'

            //    STOP 13

            //  ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //   28 RETURN

            //  END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************

            //  SUBROUTINE INTPYE (NT, S, PAR, QH, DS, XL)
            //C   Given a stress state (S) inside the PYE and a stress increment
            //C   (DS) such that (S+DS) is outside the PYE, the routine computes
            //C   the factor XL (0..1) where the stress SI=S+XL*DS, lays on the PYE.
            //C   Also, given a stress state (S) on the PYE (i.e.CC=0),
            //C   and a stress increment (DS) such that (S+DS) is outside the PYE,
            //C   but gone through unloading (i.e.crossing the inside of the PYE),
            //C   the routine computes the factor XL (0..1) where the stress SI=S+XL*DS,
            //C   is on the PYE surface.
            //C
            //C   Input: NT, S, PAR, QH, DS
            //C   Output: XL
            //C

            //  IMPLICIT REAL*8 (A-H, O-Z)

            //  DIMENSION S(1), DS(1), PAR(1), QH(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C   Initialize AA and BB

            //  AA=DABS(DS(1)*DS(1))

            //  BB=DS(1)*(S(1)-QH(NT+7))
            //C   Finalize AA and BB

            //  DO 10 I=2, NT

            //  BB=BB+(DS(I)/PAR(I+2))*((S(I)-QH(NT+I+6))/PAR(I+2))
            //   10 AA=AA+DABS(DS(I)/PAR(I+2))*DABS(DS(I)/PAR(I+2))
            //C  --------------------------------------------------------------------------------------------
            //C   Check AA

            //  IF (AA.EQ.0.D0) THEN

            //      WRITE (6,101)
            //          WRITE (11,101)
            //  101     FORMAT (1X,'SUBROUTINE INTERSECT INTPYE: VECTOR DS=0')
            //          XL=1.D0

            //      GOTO 28

            //  END IF
            //C  --------------------------------------------------------------------------------------------
            //C   Compute the Fvalue of PYE (=CC)

            //  CALL FPYE (NT, S, PAR, QH, CC)
            //C   Compute SUBSQUARE
            //C   Compute SUBSQUARE

            //  SUBSQUARE=DABS(BB*BB)-AA*CC
            //C   Compute XL

            //  XL=(-BB+DSQRT(DABS(SUBSQUARE)))/AA
            //C   Check if XL is a real number or NaN

            //  IF (DABS(XL).GE.0.D0) THEN
            //C   If F(s)=CC=0 we want the non-zero root:

            //  IF (CC.EQ.0.D0) XL=1.D0

            //  IF (XL.GT.1.D0) XL=1.D0

            //  IF (XL.LT.0.D0) XL=0.D0

            //  ELSE

            //    WRITE (6,*) 'INTPYE LINE 1161 FUCKED'

            //    STOP 14

            //  ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //   28 RETURN

            //  END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************

            //  SUBROUTINE COMPSL (NT, PAR, S, QH)
            //C   Computes the location SL of the center of the yield surface
            //C   when it is tangent to the bounding surface at the stress point (S)
            //C   Input:  NT, PAR, S, QH
            //C   Output: QH

            //  IMPLICIT REAL*8 (A-H, O-Z)

            //  DIMENSION S(1), PAR(1), QH(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C   sigma_L(1)=(1-ksi)*sigma+ksi*a

            //  QH(NT+7)=(1.D0-PAR(NT+9))*S(1)+PAR(NT+9)*QH(1)
            //C   S_L(I)=(1-ksi)*s

            //  DO 10 I=2, NT
            //   10 QH(NT+I+6)=(1.D0-PAR(NT+9))*S(I)
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************

            //  RETURN

            //  END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************

            //  SUBROUTINE COMPLDOT (NT, PAR, QH, S, DE, XLDOT)
            //C   Computes XLDOT=Lamba_dot.This routine is called only for
            //C   plastic loading (i.e.NSURF.NE.0)
            //C   Input:  NT, PAR, QH, S, DE
            //C   Output: XLDOT
            //C

            //  IMPLICIT REAL*8 (A-H, O-Z)

            //  DIMENSION QH(1), PAR(1), S(1), DE(1)

            //  DIMENSION Q(6), P(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C   Compute vector Qf at current stress state, wrt the yield surface

            //  CALL COMPQPYE (NT, PAR, QH, S, Q)
            //C   Compute array Pf

            //  CALL COMPP (NT, PAR, QH, S, P)
            //C   Compute  XN=Q*C*DE and XD=Q*C*P

            //  CALL QUAD (NT, PAR, QH, S, Q, DE, XN)

            //  CALL QUAD (NT, PAR, QH, S, Q, P, XD)
            //C   Compute elastoplastic modulus H (wrt the yield surface, i.e Hf)
            //C   The value of XD is passed to the routine because the routine
            //C   needs the sign of lambda_dot

            //  CALL COMPH (NT, PAR, QH, S, H)
            //C   Compute XD=H+Q*C*P

            //  XD=XD+H

            //  IF (DABS(XN).LE.1.D-15) THEN

            //     XLDOT=0.D0

            //  ELSE

            //     IF (XD.EQ.0.D0) XD=1.D-15

            //     XLDOT=XN/XD

            //  ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************

            //  RETURN

            //  END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************

            //  SUBROUTINE COMPYLDOT (NT, PAR, S, QH, DS, DE, YLDOT, DA)
            //C   Computes l_dot, where l_dot*beta=translation of yield surface
            //C   in the notation I am using I could be also calling it dm = l_dot
            //C Input:  NT,PAR,S,QH,DS,DA
            //C   Output: YLDOT
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION S(1), PAR(1), QH(1), DS(1), DE(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute XN, XD
            //    XN = 0.D0
            //      XD=0.D0
            //      RATIO = DA / QH(1)
            //      DO 10 I=2,NT
            //C   Compute XN = Qf:(dsigma-(a_dot/a)* sigma)
            //       XN=XN+((S(I)-QH(NT+I+6))*(DS(I)-RATIO* S(I)))/(PAR(I+2)* PAR(I+2))
            //   10 XD=XD+(S(I)-QH(NT+I+6))* S(I)/(PAR(I+2)* PAR(I+2))
            //C
            //C   Finalize XD & XN
            //      XN = XN + (S(1) - QH(NT + 7)) * (DS(1) - RATIO * S(1))
            //      XD=XD+(S(1)-QH(NT+7))*(S(1)-QH(1))
            //      XD=PAR(NT+9)*(QH(1)* QH(1))-XD
            // C   Check XD
            //      IF(XD.EQ.0.D0) XD=1.D-15
            //C  --------------------------------------------------------------------------------------------
            //C Compute l_dot=dm
            //C  --------------------------------------------------------------------------------------------
            //      YLDOT=XN/XD
            //      IF(DABS(XN).LE.1.D-15) YLDOT=0.D0
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE COMPQPYE(NT, PAR, QH, S, Q)
            //C Computes the gradient of the yield surface f.
            //C I will need to find whether I have loading or unloading.
            //C   In that sense it is obvious that the loading or unloading
            //C SHOULD be computed with respect to the PYE rather than the SSE.
            //C
            //C   Input: NT, S, PAR, QH
            //C   Output: Q
            //C

            //  IMPLICIT REAL*8 (A-H, O-Z)
            //      DIMENSION PAR(1),QH(1),S(1)
            //      DIMENSION Q(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute the volumetric part of Q
            //    Q(1)=2.D0*(S(1)-QH(NT+7))
            //C Compute the deviatoric part of Q
            //    DO 20 I=2,NT
            //   20 Q(I)=2.D0*((S(I)-QH(NT+I+6))/(PAR(I+2)* PAR(I+2)))
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE  COMPP(NT, PAR, QH, S, P)
            //C Computes the plastic potential tensor(Pf) of the PYE
            //C   I assume that it will be equal to the Q of the conjugate
            //C   point on the SSE
            //C
            //C   Input: NT,Q
            //C   Output: P
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1), QH(1), S(1)
            //      DIMENSION P(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute the Pf = Qf
            //C  --------------------------------------------------------------------------------------------
            //      CALL COMPQPYE(NT, PAR, QH, S, P)
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE  COMPPBSE(NT, PAR, QH, S, P)
            //C Computes the plastic potential tensor(Pf) of the PYE
            //C   I assume that it will be equal to the Q of the conjugate
            //C   point on the SSE
            //C
            //C   Input: NT,Q
            //C   Output: P
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1), QH(1), S(1)
            //      DIMENSION P(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute the Pf
            //C  --------------------------------------------------------------------------------------------
            //C Compute the volumetric part of P
            //    P(1)=2.D0*(S(1)-QH(1))* PAR(NT+9)
            //C Compute the deviatoric part of P
            //    DO 20 I=2,NT
            //   20 P(I)=2.D0*(S(I)/(PAR(I+2)* PAR(I+2)))* PAR(NT+9)
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE COMPDOTPDOT(NT, P, DPP)
            //C Computes the DSQRT(2/3*P':P')
            //C Input: NT,P
            //C   Output: DPP
            //C
            //      IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION P(1)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //      DPP=0.D0
            //      DO 20 I=2,NT
            //   20 DPP=DPP+DABS(P(I)* P(I))
            //      DPP=DPP/1.5D0
            //      DPP=DSQRT(DPP)
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE CONJUGATE(NT, PAR, QH, S, SM)
            //C Computes the H-conjugate stress point SH on the bounding surface,
            //C   given a stress S inside the bounding surface.The SM point lies on
            //C the SSE and on the PYE
            //C
            //C Input: NT, PAR, S, QH
            //C   Output: SH, SM
            //C
            //C The following line should be inlcuded in cyclic
            //    IMPLICIT REAL*8 (A-H, O-Z)
            //      DIMENSION PAR(1), QH(1), S(1)
            //      DIMENSION SM(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Compute the isotropic component of S_M
            //    SM(1)=QH(1)+(S(1)-QH(NT+7))/PAR(NT+9)
            //C Compute the deviatoric component of S_M
            //    DO 10 I=2,NT
            //   10 SM(I)=(S(I)-QH(NT+I+6))/PAR(NT+9)
            //C  --------------------------------------------------------------------------------------------
            //C********************************************************************************************
            //      RETURN
            //      END
            //C********************************************************************************************
            //C********************************************************************************************
            //C********************************************************************************************
            //      SUBROUTINE COMPHFACT(NT, PAR, QH, S, SM, HFACT)
            //C Computes the HFACT the factor Hf=HFACT* H(SM on the SSE)
            //C
            //C   Input: NT,PAR,S,SM,QH
            //C   Output: HFACT
            //C
            //C The following line should be inlcuded in cyclic
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1), QH(1), S(1), SM(1), Q(6), P(6)
            //C
            //C********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Initialize terms
            //    AA = 0.D0
            //      HALFTERM1=0.D0
            //      HALFTERM2 = 0.D0
            //      BB=0.D0
            //C   Start computation of AA and BB
            //      DO 31 I=2,NT
            //      HALFTERM1 = (SM(I) - S(I)) / PAR(I + 2)
            //      HALFTERM2=(S(I)-QH(I+6))/PAR(I+2)
            //      AA=AA+DABS(HALFTERM1* HALFTERM1)
            //   31 BB=BB+DABS(HALFTERM2* HALFTERM2)
            //C Finalize AA and BB
            //    AA = AA + DABS((SM(1) - S(1)) * (SM(1) - S(1)))
            //      BB=BB+DABS((S(1)-QH(7))*(S(1)-QH(7)))
            //C CHECK whether this is the first time at which I
            //C get to the plastic zone -> S==S_initial
            //    IF (DABS(AA).LE.1.D-15) THEN
            //        HFACT=0.D0
            //    ELSE
            //C Compute vector Qf at current stress state, wrt the yield surface
            //        CALL COMPQPYE (NT,PAR,QH,S,Q)
            //C Compute array Pf
            //       CALL COMPP (NT,PAR,QH,S,P)
            //C Compute XD=Q*C*P
            //       CALL QUAD (NT,PAR,QH,S,Q,P,XD)
            //         IF (DABS(BB).EQ.0.D0) BB=1.D-15
            //         HFACT=PAR(NT+10)*XD*(DABS(AA/BB)**PAR(NT+11))
            //C Check whether HFACT is a real number or NaN
            //       IF  (DABS(HFACT).GE.0.D0) THEN
            //        GOTO 34
            //        ELSE
            //        WRITE (6,*) 'COMPHFACT LINE 1443 FUCKED'
            //          STOP 15
            //       ENDIF
            //    ENDIF
            //C  --------------------------------------------------------------------------------------------
            //C  ********************************************************************************************
            //   34 RETURN
            //    END
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //      SUBROUTINE COMPH (NT,PAR,QH,S,H)
            //C Compute modulus H, wrt the yield surface. The routine is called
            //C only during plastic loading (NSURF.NE.0)
            //C
            //C Input: NT,PAR,QH,S
            //C Output: H
            //C
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1), QH(1), S(1)
            //      DIMENSION SM(6), P(6)
            //C
            //C  ********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Check if the stress point is on/inside the bounding surface
            //    IF (QH(3).EQ.1.D0) THEN
            //C  --------------------------------------------------------------------------------------------
            //C Stress point is ON the bounding surface (NSURF=1)
            //C
            //C Compute the Pf of the PYE
            //      CALL COMPP (NT,PAR,QH,S,P)
            //C Compute the DSQRT(2/3*P':P')
            //        CALL COMPDOTPDOT (NT,P,DPP)
            //C DEXPART=EXP(-(n_vp*E_vp+n_qp*E_qp))
            //        DEXPART=DEXP(-DABS(PAR(NT+5)*QH(4)+PAR(NT+6)*QH(5)))
            //C DPART1=B0-Bres+z_vp*E_vp+z_qp*E_qp
            //      DPART1=DABS(PAR(NT+3)-PAR(NT+4)+PAR(NT+7)*QH(4)+PAR(NT+8)*QH(5))
            //C DPART2=n_vp*|P|+n_qp*DSQRT(2/3*P':P')
            //        DPART2=DABS(PAR(NT+5)*DABS(P(1))+PAR(NT+6)*DPP)
            //C DPART3=z_vp*|P|+z_qp*DSQRT(2/3*P':P')
            //        DPART3=DABS(PAR(NT+7)*DABS(P(1))+PAR(NT+8)*DPP)
            //C HPART1=(1+e)/(lamda-kappa)*P*(DPART1*DEXPART+Bres)
            //        HPART1=QH(2)/(PAR(2)-PAR(3))*P(1)*(DPART1*DEXPART+PAR(NT+4))
            //C HPART2=DPART3*DEXPART
            //      HPART2=DABS(DPART3*DEXPART)
            //C HPART3=-DPART1*DPART2*DEXPART
            //      HPART3=-DABS(DPART1*DPART2*DEXPART)
            //C Compute the plastic Modulus H
            //      H=2.D0*S(1)*QH(6)*PAR(NT+9)*(HPART1+HPART2+HPART3)
            //C  --------------------------------------------------------------------------------------------
            //      ELSE
            //C  --------------------------------------------------------------------------------------------
            //C Stress point is inside the bounding surface (NSURF=-1)
            //C RECompute the conjugate point SM according to Kavvadas et. al.
            //        CALL CONJUGATE (NT,PAR,QH,S,SM)
            //C Compute HFACT -> Hf=HFACT*H(SM on the SSE)
            //        CALL COMPHFACT (NT,PAR,QH,S,SM,HFACT)
            //C Compute the Pf of the conjugate point on SSE
            //      CALL COMPPBSE (NT,PAR,QH,SM,P)
            //C Compute the DSQRT(2/3*P':P')
            //        CALL COMPDOTPDOT (NT,P,DPP)
            //C DEXPART=EXP(-(n_vp*E_vp+n_qp*E_qp))
            //        DEXPART=DEXP(-(PAR(NT+5)*QH(4)+PAR(NT+6)*QH(5)))
            //C DPART1=B0-Bres+z_vp*E_vp+z_qp*E_qp
            //      DPART1=PAR(NT+3)-PAR(NT+4)+PAR(NT+7)*QH(4)+PAR(NT+8)*QH(5)
            //C DPART2=n_vp*|P|+n_qp*DSQRT(2/3*P':P')
            //        DPART2=PAR(NT+5)*DABS(P(1))+PAR(NT+6)*DPP
            //C DPART3=z_vp*|P|+z_qp*DSQRT(2/3*P':P')
            //        DPART3=PAR(NT+7)*DABS(P(1))+PAR(NT+8)*DPP
            //C HPART1=(1+e)/(lamda-kappa)*P*(DPART1*DEXPART+Bres)
            //        HPART1=QH(2)/(PAR(2)-PAR(3))*P(1)*(DPART1*DEXPART+PAR(NT+4))
            //C HPART2=DPART3*DEXPART
            //      HPART2=DPART3*DEXPART
            //C HPART3=-DPART1*DPART2*DEXPART
            //      HPART3=-DPART1*DPART2*DEXPART
            //C Compute the plastic Modulus H
            //      Hf1=2.D0*SM(1)*QH(6)*PAR(NT+9)*(HPART1+HPART2+HPART3)
            //        Hf2=HFACT
            //      H=Hf1+Hf2
            //      QH(2*NT+7)=Hf1
            //      QH(2*NT+8)=Hf2
            //C H=Hf1+HFACT
            //C WRITE (11,*) H, Hf1, HFACT
            //C WRITE (11,*) H
            //C  --------------------------------------------------------------------------------------------
            //      END IF
            //C  --------------------------------------------------------------------------------------------
            //C  ********************************************************************************************
            //      RETURN
            //    END
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //      SUBROUTINE CSTIF (ISTIF,NT,QH,S,PAR,C)
            //C Computes the elastic OR elastoplastic stiffness C(NT,NT)
            //C at the current state (S,QH)
            //C          if ISTIF=0: compute elastic stiffness
            //C          if ISTIF=1: compute elasto-plastic stiffness
            //C Two types of elasticity are implemented:
            //C IMODE=1: poroelasticity (Cam-Clay)
            //C Input:  ISTIF,NT,QH,S,PAR
            //C Output: C
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION QH(1),S(1),PAR(1),C(NT,NT)
            //      DIMENSION Q(6),P(6),DD(6),EE(6)
            //C
            //C  ********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Zero stiffness C:
            //      DO 1 I=1,NT
            //    DO 1 J=1,NT
            //    1 C(I,J)=0.D0
            //C Compute parameters
            //C Poro-elasticity (Cam-Clay): COMPUTE K, 2G
            //       XK=QH(2)*S(1)/PAR(3)
            //         TWOG=PAR(1)*XK
            //C  --------------------------------------------------------------------------------------------
            //C Compute elastic stiffness
            //    DO 60 I=2,NT
            //   60 C(I,I)=TWOG
            //C Poro-elasticity
            //    C(1,1)=XK
            //C  ********************************************************************************************
            //C  *****************************************  PROSOXH  ****************************************
            //C  ********************************************************************************************
            //C  *****************************************  RETURN  *****************************************
            //C If only elastic stiffness is required: finished
            //    IF (ISTIF.EQ.0) RETURN
            //C  --------------------------------------------------------------------------------------------
            //C Compute elastoplastic stiffness:
            //C Compute vectors P,Q wrt the yield surface (assumes P=Q)
            //C
            //C Compute vector Qf at current stress state, wrt the yield surface
            //    CALL COMPQPYE (NT,PAR,QH,S,Q)
            //C Compute array Pf
            //    CALL COMPP (NT,PAR,QH,S,P)
            //C Compute DD=C*P
            //    CALL COMPDS (NT,PAR,QH,S,P,DD)
            //C Compute temporary OMEGA=Q*C*P
            //    OMEGA=0.D0
            //    DO 115 I=1,NT
            //  115 OMEGA=OMEGA+Q(I)*DD(I)
            //C Compute elastoplastic modulus H (wrt the yield surface)
            //      CALL COMPH (NT,PAR,QH,S,H)
            //C Compute final OMEGA=H+Q*C*P
            //    OMEGA=H+OMEGA
            //    IF (OMEGA.EQ.0.D0)  THEN
            //      WRITE (6,*) 'SUBROUTINE CSTIF FAILED - OMEGA = 0.D0'
            //        STOP 6
            //      ENDIF
            //C Compute EE=Q*C
            //    CALL COMPDS (NT,PAR,QH,S,Q,EE)
            //C Compute C = Celastic - (C*P)(Q*C)/OMEGA
            //    DO 120 I=1,NT
            //    DO 120 J=1,NT
            //  120 C(I,J)=C(I,J)-DD(I)*EE(J)/OMEGA
            //C  --------------------------------------------------------------------------------------------
            //C  ********************************************************************************************
            //      RETURN
            //    END
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //      SUBROUTINE COMPDS (NT,PAR,QH,S,DE,DS)
            //C
            //C Computes the stress increment DS=C*DE for an elastic strain increment
            //C DE (both in the transformed space). C=elasticity matrix.
            //C Uses poro-elasticity (if IMODE=1) [or hyper-elasticity(if IMODE = 2)-Not for CAM-CLAY]
            //    C   for the calculation of the elasticity matrix.
            //C INPUT:    NT,PAR,QH,S,DE
            //C OUTPUT:   DS(NT)
            //C
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1),QH(1),S(1),DE(1),DS(1)
            //C
            //C  ********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C Poro-elasticity
            //C Compute elastic parameters K, 2G
            //    XK=QH(2)*S(1)/PAR(3)
            //      TWOG=PAR(1)*XK
            //C Compute DS
            //    DS(1)=XK*DE(1)
            //      DO 20 I=2,NT
            //   20 DS(I)=TWOG*DE(I)
            //C  --------------------------------------------------------------------------------------------
            //C  ********************************************************************************************
            //      RETURN
            //    END
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //      SUBROUTINE QUAD (NT,PAR,QH,S,Q,P,VALUE)
            //C
            //C Forms the quadratic form Q*C*P, where Q,P are vectors and C is the
            //C elasticity matrix. The elasticity of the material is either
            //C poro-elastic (IMODE=1) or hyperelastic (IMODE=2). Returns the
            //C value (VALUE) of the quadratic form.
            //C INPUT:  NT,PAR,QH,S,Q,P
            //C OUTPUT: VALUE
            //C
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION PAR(1),QH(1),S(1),Q(1),P(1)
            //      DIMENSION DD(6)
            //C
            //C  ********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //      CALL COMPDS (NT,PAR,QH,S,P,DD)
            //      VALUE=0.D0
            //    DO 10 I=1,NT
            //   10 VALUE=VALUE+Q(I)*DD(I)
            //C  --------------------------------------------------------------------------------------------
            //C  ********************************************************************************************
            //      RETURN
            //    END
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //      SUBROUTINE TRNSFRM (IMODE,NT,S,SR,XJ,XJR)
            //C
            //C Transforms stress, strain and Jacobian matrix between standard
            //C and transformed spaces.
            //C Standard:     X, Y, Z, XY, XZ, YZ (or R, Z, TH, RZ, RTH, ZTH)
            //C Transformed:  Greek, 1, 2, 3, 4, 5
            //C NT= dimension of space ( 3<= NT <= 6 )
            //C NT=3  : only principal stresses/strains
            //C      4  : plane strain / axisymmetric
            //C      6  : general 3-D
            //C
            //C IMODE=  1: Stress,    Standard     ->  Transformed : S(tr)=Ts    * S
            //C           2: Stress,    Transformed  ->  Standard    : S    =Ts(-1)* S(tr)
            //C           3: Strain,    Standard     ->  Transformed : E(tr)=Te    * E
            //C           4: Strain,    Transformed  ->  Standard    : E    =Te(-1)* E(tr)
            //C           5: Jacobian,  Standard     ->  Transformed (NOT IMPLEMENTED)
            //C           6: Jacobian,  Transformed  ->  Standard    : J = Te(T)*J(tr)*Te
            //C Note that: Ts=Te(-T), Ts(-1)=Te(T)
            //C
            //C S(NT)     = Input stress or strain vector.
            //C SR(NT)    = Output stress or strain vector
            //C XJ(NT,NT)  = Input Jacobian matrix
            //C XJR(NT,NT) = Output Jacobian matrix
            //C
            //    IMPLICIT REAL*8 (A-H,O-Z)
            //      DIMENSION S(1), SR(1), XJ(NT,1), XJR(NT,1)
            //      DIMENSION A(3,3), B(3,3)
            //C
            //C  ********************************************************************************************
            //C  --------------------------------------------------------------------------------------------
            //C X1=1/SQRT(6), X2=1/SQRT(2), X3=1/3, X4=SQRT(2)
            //      X1=0.408248290464D0
            //      X2=0.707106781188D0
            //      X3=0.333333333333D0
            //      X4=1.41421356237D0
            //C  --------------------------------------------------------------------------------------------
            //C Set-up array A(3,3)
            //C
            //    A(2,1)=-X1
            //    A(3,1)=-X2
            //    A(2,2)=2.D0*X1
            //    A(3,2)=0.D0
            //    A(2,3)=-X1
            //    A(3,3)=X2
            //C
            //    GOTO (2,5,5,2,1,5  ), IMODE
            //1     WRITE (6,98)
            //98    FORMAT (1X,'NOT IMPLEMENTED TRANSFORMATION CALL')
            //      STOP 8
            //2     CONTINUE
            //C A(3,3) stores A(-T),  XL=sqrt(2)
            //      XL=X4
            //    DO 3 I=1,3
            //3     A(1,I)=X3
            //    GOTO 9
            //5     CONTINUE
            //C A(3,3) stores A,  XL=1/sqrt(2)
            //      XL=X2
            //    DO 6 I=1,3
            //6     A(1,I)=1.D0
            //9     CONTINUE
            //C  --------------------------------------------------------------------------------------------
            //      GOTO (10,20,10,20,11,100  ), IMODE
            //11    RETURN
            //C Stress or Strain:    Standard     ->  Transformed
            //10    CONTINUE
            //    DO 12 I=1,3
            //      SR(I)=0.D0
            //    DO 12 J=1,3
            //12    SR(I)=SR(I)+A(I,J)*S(J)
            //      IF (NT.EQ.3) RETURN
            //    GOTO 50
            //C Stress or Strain:    Transformed  ->  Standard
            //20    CONTINUE
            //    DO 22 I=1,3
            //      SR(I)=0.D0
            //    DO 22 J=1,3
            //22    SR(I)=SR(I)+A(J,I)*S(J)
            //      IF (NT.EQ.3) RETURN
            //50    CONTINUE
            //    DO 53 I=4,NT
            //53    SR(I)=XL*S(I)
            //      RETURN
            //C Jacobian:  Transformed  ->  Standard
            //100   CONTINUE
            //C   (a)  Form upper left section:
            //      DO 110 I=1,3
            //      DO 110 L=1,3
            //      B(I,L)=0.D0
            //    DO 110 K=1,3
            //110   B(I,L)=B(I,L)+A(K,I)*XJ(K,L)
            //      DO 120 I=1,3
            //      DO 120 J=1,3
            //      XJR(I,J)=0.D0
            //    DO 120 L=1,3
            //120   XJR(I,J)=XJR(I,J)+B(I,L)*A(L,J)
            //      IF (NT.EQ.3) RETURN
            //C   (b)  Form lower left section:
            //      DO 130 I=4,NT
            //    DO 130 J=1,3
            //      DUM=0.D0
            //    DO 125 K=1,3
            //125   DUM=DUM+XJ(I,K)*A(K,J)
            //130   XJR(I,J)=DUM*XL
            //C   (c)  Form upper right section:
            //      DO 140 I=1,3
            //      DO 140 J=4,NT
            //    DUM=0.D0
            //    DO 135 K=1,3
            //135   DUM=DUM+A(K,I)*XJ(K,J)
            //140   XJR(I,J)=DUM*XL
            //C   (d)  Form lower right section:
            //      DUM=XL*XL
            //    DO 150 I=4,NT
            //    DO 150 J=4,NT
            //150   XJR(I,J)=DUM*XJ(I,J)
            //C  --------------------------------------------------------------------------------------------
            //C  ********************************************************************************************
            //      RETURN
            //    END
            //C  ********************************************************************************************
            //C  ********************************************************************************************
            //C  ********************************************************************************************

            #endregion
            #region TransformationForStresses
        public double[] TransformToTransformedStresses( double[] Stresses)
        {
            //Transform at first to the Soil Mechanics sign convention
            for (int i=0;i<6;i++)
            {
                Stresses[i] = -Stresses[i];
            }
            var help = new double[6];
            help[0] = (Stresses[0] + Stresses[1] + Stresses[2]) / 3;
            help[1] = (2 * Stresses[1] - Stresses[0] - Stresses[2]) / Math.Sqrt(6);
            help[2] = (Stresses[2] - Stresses[0]) / Math.Sqrt(2);
            help[3] = (Stresses[3]) / Math.Sqrt(2);
            help[4] = (Stresses[4]) / Math.Sqrt(2);
            help[5] = (Stresses[5]) / Math.Sqrt(2);
            Stresses = help;
            return Stresses;
        }
        public double[] TransformToStandardStresses(double[] Stresses)
        {
            var help = new double[6];
            help[0] = Stresses[0] - Stresses[1] / Math.Sqrt(6) - Stresses[2] / Math.Sqrt(2);
            help[1] = Stresses[0] + 2 * Stresses[1] / Math.Sqrt(6);
            help[2] = Stresses[0] - Stresses[1] / Math.Sqrt(6) + Stresses[2] / Math.Sqrt(2);
            help[3] = (Stresses[3]) * Math.Sqrt(2);
            help[4] = (Stresses[4]) * Math.Sqrt(2);
            help[5] = (Stresses[5]) * Math.Sqrt(2);
            Stresses = help;
            //Transform at the end to the engineering convention
            for (int i = 0; i < 6; i++)
            {
                Stresses[i] = -Stresses[i];
            }
            return Stresses;
        }
        #endregion
        #region TransformationForStrains
        //public double[] TransformToTransformedStrains(out double[] Strains)
        //{
        //    Transform at first to the Soil Mechanics sign convention
        //    for (int i=0;i<6;i++)
        //    {
        //       Strains[i] = -Strains[i];
        //    }
        //    var help = new double[6];
        //    help[0] = (Strains[0] + Strains[1] + Strains[2]);
        //    help[1] = (2 * Strains[1] - Strains[0] - Strains[2]) / Math.Sqrt(6);
        //    help[2] = (Strains[2] - Strains[0]) / Math.Sqrt(2);
        //    help[3] = (Strains[3]) / Math.Sqrt(2);
        //    help[4] = (Strains[4]) / Math.Sqrt(2);
        //    help[5] = (Strains[5]) / Math.Sqrt(2);
        //    Strains = help;
        //    return Strains;
        //}
        //public double[] TransformToStandardStrains(out double[] Strains)
        //{
        //    var help = new double[6];
        //    help[0] = Strains[0]/3 - Strains[1] / Math.Sqrt(6) - Strains[2] / Math.Sqrt(2);
        //    help[1] = Strains[0]/3 + 2 * Strains[1] / Math.Sqrt(6);
        //    help[2] = Strains[0]/3 - Strains[1] / Math.Sqrt(6) + Strains[2] / Math.Sqrt(2);
        //    help[3] = (Strains[3]) * Math.Sqrt(2);
        //    help[4] = (Strains[4]) * Math.Sqrt(2);
        //    help[5] = (Strains[5]) * Math.Sqrt(2);
        //    Strains = help;
       //    Transform at the end to the engineering convention
       //     for (int i = 0; i< 6; i++)
       //     {
       //        Strains[i] = -Strains[i];
       //     }
    //    return Strains;
    //}
    #endregion
    #region TransformationForJacobianMatrix 
    //We consider the case of associative plasticity. In other words the plastic potential is the Yield Surface. Therefore the Consistent Constitutive Matrix is symmetric. 
    //The Matrix Mutiplication is done directly from Matlab
    public IMatrix2D<double> TransformationForJacobianMatrix (IMatrix2D<double> ConstitutiveMatrix)
        {
            var m1 = new double[21];
            var s1 = new double[25];
            //Helping assigniments
            ConstitutiveMatrix[0, 0] = m1[0];
            ConstitutiveMatrix[0, 1] = m1[1];
            ConstitutiveMatrix[0, 2] = m1[2];
            ConstitutiveMatrix[0, 3] = m1[3];
            ConstitutiveMatrix[0, 4] = m1[4];
            ConstitutiveMatrix[0, 5] = m1[5];
            ConstitutiveMatrix[1, 1] = m1[6];
            ConstitutiveMatrix[1, 2] = m1[7];
            ConstitutiveMatrix[1, 3] = m1[8];
            ConstitutiveMatrix[1, 4] = m1[9];
            ConstitutiveMatrix[1, 5] = m1[10];
            ConstitutiveMatrix[2, 2] = m1[11];
            ConstitutiveMatrix[2, 3] = m1[12];
            ConstitutiveMatrix[2, 4] = m1[13];
            ConstitutiveMatrix[2, 5] = m1[14];
            ConstitutiveMatrix[3, 3] = m1[15];
            ConstitutiveMatrix[3, 4] = m1[16];
            ConstitutiveMatrix[3, 5] = m1[17];
            ConstitutiveMatrix[4, 4] = m1[18];
            ConstitutiveMatrix[4, 5] = m1[19];
            ConstitutiveMatrix[5, 5] = m1[20];
            s1[0] =m1[1]/Math.Sqrt(6);
            s1[1] =m1[2]/Math.Sqrt(2);
            s1[2] =(m1[2]+Math.Sqrt(6)*m1[7]/3)/Math.Sqrt(2);
            s1[3] = m1[1] * Math.Sqrt(6)/3 ;
            s1[4] = m1[5] / Math.Sqrt(2); 
            s1[5] = m1[4] / Math.Sqrt(2);
            s1[6] = m1[3] / Math.Sqrt(2);
            s1[7] = Math.Sqrt(6) * (m1[1]+m1[6]*Math.Sqrt(6)/3);
            s1[8] = Math.Sqrt(12) * m1[10];
            s1[9] = Math.Sqrt(12) * m1[9];
            s1[10] =Math.Sqrt(12) * m1[8];
            s1[11] =m1[10]/ Math.Sqrt(6);
            s1[12] =m1[14]/ Math.Sqrt(2);
            s1[13] =m1[9]/ Math.Sqrt(6);
            s1[14] =m1[13]/ Math.Sqrt(2);
            s1[15] =m1[8]/ Math.Sqrt(6);
            s1[16] =m1[12]/ Math.Sqrt(2);
            s1[21] = m1[7]/ Math.Sqrt(6);
            s1[22] =m1[11]/ Math.Sqrt(2);
            s1[23] = m1[6] / Math.Sqrt(6);
            s1[24] = m1[7] / Math.Sqrt(2);
            s1[17] = (s1[21]+s1[22]-m1[2]) / Math.Sqrt(2);
            s1[18] = (-s1[21] + s1[22] + m1[2]) / Math.Sqrt(2);
            s1[19] = Math.Sqrt(6) * (s1[24]+s1[23]-m1[1]);
            s1[20] = Math.Sqrt(6) * (s1[24] - s1[23] + m1[1]);
            //Transform
            ConstitutiveMatrix[0, 0] = m1[0]-s1[0]-s1[1]+s1[19]/6+s1[17];
            ConstitutiveMatrix[0, 1] = m1[0] - s1[0] - s1[1] - s1[19] / 3 ;
            ConstitutiveMatrix[0, 2] = m1[0] - s1[0] - s1[1] + s1[19] / 6 - s1[17];
            ConstitutiveMatrix[0, 3] = -(s1[15]+s1[16]-m1[3])/Math.Sqrt(2);
            ConstitutiveMatrix[0, 4] = -(s1[13] + s1[14] - m1[4]) / Math.Sqrt(2);
            ConstitutiveMatrix[0, 5] = -(s1[11] + s1[12] - m1[5]) / Math.Sqrt(2);
            ConstitutiveMatrix[1, 1] = m1[0] + s1[3] + s1[7] / 3;
            ConstitutiveMatrix[1, 2] = m1[0] + s1[3] - s1[7] / 6 + s1[2];
            ConstitutiveMatrix[1, 3] = s1[6] + s1[10] / 6;
            ConstitutiveMatrix[1, 4] = s1[5] + s1[9] / 6;
            ConstitutiveMatrix[1, 5] = s1[4] + s1[8] / 6;
            ConstitutiveMatrix[2, 2] = m1[0] + s1[1] - s1[0] - s1[20] / 6 + s1[18];
            ConstitutiveMatrix[2, 3] = (s1[16] - s1[15] + m1[3]) / Math.Sqrt(2);
            ConstitutiveMatrix[2, 4] = (s1[14] - s1[13] + m1[4]) / Math.Sqrt(2);
            ConstitutiveMatrix[2, 5] = (s1[12] - s1[11] + m1[5]) / Math.Sqrt(2);
            ConstitutiveMatrix[3, 3] = m1[15]/2;
            ConstitutiveMatrix[3, 4] = m1[16]/2;
            ConstitutiveMatrix[3, 5] = m1[17]/2;
            ConstitutiveMatrix[4, 4] = m1[18]/2;
            ConstitutiveMatrix[4, 5] = m1[19]/2;
            ConstitutiveMatrix[5, 5] = m1[20]/2;
            //Apply Symmetry      
            ConstitutiveMatrix[1, 0] = ConstitutiveMatrix[0, 1];
            ConstitutiveMatrix[2, 0] = ConstitutiveMatrix[0, 2];
            ConstitutiveMatrix[3, 0] = ConstitutiveMatrix[0, 3];
            ConstitutiveMatrix[4, 0] = ConstitutiveMatrix[0, 4];
            ConstitutiveMatrix[5, 0] = ConstitutiveMatrix[0, 5];
            ConstitutiveMatrix[2, 1] = ConstitutiveMatrix[1, 2];
            ConstitutiveMatrix[3, 1] = ConstitutiveMatrix[1, 3];
            ConstitutiveMatrix[4, 1] = ConstitutiveMatrix[1, 4];
            ConstitutiveMatrix[5, 1] = ConstitutiveMatrix[1, 5];
            ConstitutiveMatrix[3, 2] = ConstitutiveMatrix[2, 3];
            ConstitutiveMatrix[4, 2] = ConstitutiveMatrix[2, 4];
            ConstitutiveMatrix[5, 2] = ConstitutiveMatrix[2, 5];
            ConstitutiveMatrix[4, 3] = ConstitutiveMatrix[3, 4];
            ConstitutiveMatrix[5, 3] = ConstitutiveMatrix[3, 5];
            ConstitutiveMatrix[5, 4] = ConstitutiveMatrix[4, 5];
            this.ConstitutiveMatrix = ConstitutiveMatrix;
            return ConstitutiveMatrix;
        }
        #endregion
        #region LoadSubroutines
        public void DecideLoad(double[] PAR,double [] de, double[] QH, double[] Stresses, double[] elasticConstitutiveMatrix, IMatrix2D<double> ConstitutiveMatrix)
        {
            Stresses = TransformToTransformedStresses(Stresses);
            if (QH[3] == 0)
            {
                Stresses=unload( PAR, de,  QH,  Stresses);
            }
            else
            {
                var help = 0.0;
                help = Numer(PAR, de, QH, Stresses);
                if (help < 0)
                {
                    Stresses=unload( PAR, de,  QH,  Stresses);
                   }
                else if (help > 0)
                {
                    Stresses=load(PAR, de, QH, Stresses);
                    }
                else if (help == 0)
                {
                    QH=loadneutral(PAR, de, QH, Stresses);
                   }
            }
            ConstitutiveMatrix = cstif(PAR, QH, Stresses);
            ConstitutiveMatrix = TransformationForJacobianMatrix(ConstitutiveMatrix);
            this.Stresses = Stresses;
            Stresses = TransformToStandardStresses(Stresses);
        }
        public double [] load(double[] PAR,double [] de, double[] QH, double[] Stresses)
        {
            if (QH[3] == -1)
            {
                Stresses=yload( PAR, de, QH, Stresses);
            }
            else
            {
                Stresses=bload( PAR, de, QH, Stresses);
            }
            return Stresses;
        }
        public double[] bload( double[] PAR,double [] de, double[] QH, double[] Stresses)
        {
            var Dpp = 0.0;
            var XLDOT = 0.0;
            var Devp = 0.0;
            var Deqp = 0.0;
            var Devptrue = 0.0;
            var Compevp = 0.0;
            var Compeqp = 0.0;
            var Dexpart = 0.0;
            var Dpart1 = 0.0;
            var DD = new double[6];
            var ds = new double[6];
            var Alpha = 0.0;
            XLDOT = compldot(PAR, de, QH, Stresses);
            Qgrad = compp(PAR, QH, Stresses);
            Dpp = compdotpdot(Qgrad);
            Devp = Math.Abs(XLDOT * Qgrad[0]);
            Deqp = Math.Abs(XLDOT * Dpp);
            Devptrue = Math.Abs(XLDOT) * Qgrad[0];
            Compevp = QH[4] + Devp;
            Compeqp = QH[5] + Deqp;
            Dexpart = Math.Exp(-Math.Abs(PAR[11] * Compevp + PAR[12] * Compeqp));
            Dpart1 = Math.Abs(PAR[9] - PAR[10] + PAR[13] * Compevp + PAR[14] * Compeqp);
            QH[4] = Compevp;
            QH[5] = Compeqp;
            QH[6] = QH[6] * (1 + QH[2] * Devptrue / (PAR[2] - PAR[3]));
            QH[1] = QH[6] * (Dpart1 * Dexpart + PAR[10]);
            for (int i = 0; i < 6; i++)
            {
                DD[i] = de[i] - Math.Abs(XLDOT) * Qgrad[i];
            }
            ds = compds(PAR, DD, QH, Stresses);
            for (int i = 0; i < 6; i++)
            {
                Stresses[i] = Stresses[i] + ds[i];
            }
            Alpha = QH[1];
            QH[2] = QH[2] * (1 - de[0]);
            Stresses=adjbse(PAR, Alpha, Stresses);
            compsl(PAR, QH, Stresses);
            QH[3] = 1;
            this.QH = QH;
            return Stresses;
        }
        public double[] yload( double[] PAR,double [] de,  double[] QH,  double[] Stresses)
        {
            var Dpp = 0.0;
            var XLDOT = 0.0;
            var Devp = 0.0;
            var Deqp = 0.0;
            var Devptrue = 0.0;
            var Compevp = 0.0;
            var Compeqp = 0.0;
            var Dexpart = 0.0;
            var Dpart1 = 0.0;
            var QH6 = 0.0;
            var Alpha = 0.0;
            var DA = 0.0;
            var Ratio = 0.0;
            var YPLDOT = 0.0;
            var YLDOT = 0.0;
            var Fcheck = 0.0;
            var RR = 0.0;
            var B = new double[6];
            var ds = new double[6];
            var DD = new double[6];
            var DD1 = new double[6];
            var Stressescheck = new double[6];
            XLDOT = compldot(PAR, de, QH, Stresses);
            Qgrad = compp(PAR, QH, Stresses);
            Dpp = compdotpdot(Qgrad);
            for (int i = 0; i < 6; i++)
            {
                DD[i] = de[i] - Math.Abs(XLDOT) * Qgrad[i];
            }
            ds = compds(PAR, DD, QH, Stresses);
            B[0] = (Stresses[0] - QH[13]) / PAR[15] - (Stresses[0] - QH[1]);
            for (int i = 1; i < 6; i++)
            {
                B[i] = (Stresses[i] - QH[i + 13]) / PAR[15] - Stresses[i];
            }
            Devp = Math.Abs(XLDOT * Qgrad[0]);
            Deqp = Math.Abs(XLDOT * Dpp);
            Devptrue = Math.Abs(XLDOT) * Qgrad[0];
            Compevp = QH[4] + Devp;
            Compeqp = QH[5] + Deqp;
            Dexpart = Math.Exp(-Math.Abs(PAR[11] * Compevp + PAR[12] * Compeqp));
            Dpart1 = Math.Abs(PAR[9] - PAR[10] + PAR[13] * Compevp + PAR[14] * Compeqp);
            QH6 = QH[6] * (1 + QH[2] * Devptrue / (PAR[2] - PAR[3]));
            Alpha = QH6 * (Dpart1 * Dexpart + PAR[10]);
            DA = Alpha - QH[1];
            Ratio = DA / QH[1];
            YPLDOT = 1 + Ratio;
            for (int i = 0; i < 6; i++)
            {
                Stressescheck[i] = Stresses[i] + ds[i];
            }
            Fcheck = FBSE(PAR, Alpha, Stressescheck);
            YLDOT = compyldot(PAR, QH, Stresses, ds, de, DA);
            YPLDOT = 1 + Ratio;
            if (Fcheck<0)
            {
                QH[1] = Alpha;
                for (int i = 0; i < 6; i++)
                {
                    QH[13 + i] = YPLDOT * QH[13 + i] + YLDOT * B[i];          
                }
                for (int i = 0; i < 6; i++)
                {
                    Stresses[i] = Stresses[i] + ds[i];
                }
                QH[4] = Compevp;
                QH[5] = Compeqp;
                QH[6] = QH6;
                QH[2] = QH[2] * (1 - de[0]);
                Stresses=adjpye(PAR,QH,Stresses);
                QH[3] = -1;
            }
            else
            {
                RR = IntBse(PAR, QH, Stresses, ds);
                for (int i = 0; i < 6; i++)
                {
                    DD[i] = (1- RR)*de[i];
                    DD1[i] = (RR) * de[i];
                }
                XLDOT = compldot(PAR, DD1, QH, Stresses);
                Qgrad = compp(PAR, QH, Stresses);
                Dpp = compdotpdot(Qgrad);
                Devp = Math.Abs(XLDOT * Qgrad[0]);
                Deqp = Math.Abs(XLDOT * Dpp);
                Devptrue = Math.Abs(XLDOT) * Qgrad[0];
                Compevp = QH[4] + Devp;
                Compeqp = QH[5] + Deqp;
                Dexpart = Math.Exp(-Math.Abs(PAR[11] * Compevp + PAR[12] * Compeqp));
                Dpart1 = Math.Abs(PAR[9] - PAR[10] + PAR[13] * Compevp + PAR[14] * Compeqp);
                QH[6] = QH[6] * (1 + QH[2] * Devptrue / (PAR[2] - PAR[3]));
                QH[1] = QH[6] * (Dpart1*Dexpart+PAR[10]);
                QH[4] = Compevp;
                QH[5] = Compeqp;
                QH[2] = QH[2] * (1 - DD1[0]);
                for (int i = 0; i < 6; i++)
                {
                    Stresses[i] = Stresses[i] + RR*ds[i];
                }
                Stresses=adjbse(PAR, QH[1], Stresses);
                compsl(PAR, QH, Stresses);
                QH[3] = 1;
                if (RR<1)
                {
                    bload(PAR,DD, QH, Stresses);
                }
            }
            this.QH = QH;
            return Stresses;
        }
        public double [] unload( double[] PAR, double[] de,  double[] QH,  double[] Stresses)
        {
            var ds = new double[6];
            var DV = 0.0;
            var Stressescheck = new double[6];
            var Alpha = 0.0;
            var Fcheck = 0.0;
            var Fcheck2 = 0.0;
            var XLamda = 0.0;
            var DD1 = new double[6];
            ds=compds (PAR, de, QH, Stresses);
            DV = -de[0] * QH[2];
            for (int i = 0; i < 6; i++)
            {
                Stressescheck[i] = Stresses[i] + ds[i];
            }
            Alpha = QH[1];
            Fcheck = FPYE(PAR, QH, Stressescheck);
            if (Fcheck<0)
            {
                Stresses = Stressescheck;
                QH[2] = QH[2] + DV;
                QH[3] = 0;
            }
            else
            {
                XLamda=IntPye(PAR, QH, ds, Stresses);
                for (int i = 0; i < 6; i++)
                {
                    Stresses[i] = Stresses[i] + XLamda*ds[i];
                    DD1[i] = (1 - XLamda) * de[i];
                }
                QH[2] = QH[2] + XLamda * DV;
                Stresses=adjpye(PAR, QH, Stresses);
                Fcheck2 = FBSE(PAR,Alpha,Stresses);
                if (Fcheck2<0)
                {
                    QH[3] = -1;
                }
                else
                {
                    Stresses=adjbse(PAR, Alpha, Stresses);
                    compsl(PAR, QH, Stresses);
                    QH[3] = 1;
                }
                for (int i = 0; i < 6; i++)
                {
                    QH[7 + i] = Stresses[i];
                }
                load(PAR, DD1, QH, Stresses);
            }
            this.QH = QH;
            return Stresses;
        }
        public double [] loadneutral(double[] PAR,double[] de, double[] QH, double[] Stresses)
        {
            QH[2] = QH[2] * (1 - de[0]);
            return QH;
        }
        #endregion
        #region Numer/FB/FP/INT/ADJSubroutines
        public double Numer (double[] PAR, double [] de,double [] QH,double[] Stresses)
        {
            var Xnum = 0.0;
            Qgrad = compp(PAR, QH, Stresses);
            Xnum = Quad(PAR,de,QH,Stresses,Qgrad);
            return Xnum;
        }
        public double FBSE(double[] PAR,double Alpha, double[] Stresses)
        {
            var F = 0.0;
            F = Math.Pow((Alpha - Stresses[0]), 2) - Math.Pow(Alpha, 2);
            for (int i = 1; i < 6; i++)
            {
                F = F + Math.Pow((Stresses[i] / PAR[i + 3]), 2);
            }
            return F;
        }
        public double FPYE(double[] PAR, double [] QH, double[] Stresses)
        {
            var F = 0.0;
            F = Math.Pow((QH[13] - Stresses[0]), 2) - Math.Pow((PAR[15]*QH[1]), 2);
            for (int i = 1; i < 6; i++)
            {
                F = F + Math.Pow(((Stresses[i]-QH[13+i]) / PAR[i + 3]), 2);
            }
            return F;
        }
        public double[] adjbse(double [] PAR,double Alpha, double [] Stresses)
        {
            var F1 = 0.0;
            var Alpha1 = 0.0;
            var XLamda = 0.0;
            F1 = FBSE(PAR, Alpha, Stresses);
            Alpha1 = Math.Pow(Alpha, 2);
            XLamda = Math.Sqrt(Math.Abs(Alpha1 / (F1 + Alpha1)));
            Stresses[0] = Alpha + XLamda * (Stresses[0] - Alpha);
            for (int i = 1; i < 6; i++)
            {
                Stresses[i] = XLamda * Stresses[i];
            }
            return Stresses;
        }
        public double[] adjpye(double [] PAR,double [] QH, double [] Stresses)
        {
            var F1 = 0.0;
            var Alpha1 = 0.0;
            var XLamda = 0.0;
            F1 = FPYE(PAR, QH, Stresses);
            Alpha1 = Math.Pow((QH[1] * PAR[15]), 2);
            XLamda = Math.Sqrt(Math.Abs(Alpha1 / (F1 + Alpha1)));
            for (int i = 0; i < 6; i++)
            {
                Stresses[i] = XLamda * Stresses[i] + (1 - XLamda) * QH[13 + i];
            }
            return Stresses;
        }
        public double IntBse(double [] PAR,double [] QH, double [] Stresses,double [] Ds)
        {
            var AA = 0.0;
            var BB = 0.0;
            var CC = 0.0;
            var dummy = 0.0;
            var XLamda = 0.0;
            AA = Math.Abs(Math.Pow(Ds[0], 2));
            BB = Ds[0] * (Stresses[0] - QH[1]);
            for (int i = 1; i < 6; i++)
            {
                BB = BB + (Ds[i] * Stresses[i]) / Math.Pow(PAR[3 + i], 2);
                AA = AA + Math.Pow((Ds[i] / PAR[3 + i]), 2);
            }
            if (AA == 0)
            {
                XLamda = 1;
            }
            CC = FBSE(PAR, QH[1], Stresses);
            dummy = Math.Pow(BB, 2) - AA * CC;
            XLamda = (-BB + Math.Sqrt(dummy)) / (AA);
            if (CC == 0)
            {
                XLamda = 1;
            }
            if (XLamda > 1)
            {
                XLamda = 1;
            }
            if (XLamda < 0)
            {
                XLamda = 0;
            }
            return XLamda;
        }
        public double IntPye(double [] PAR,double[] QH, double [] Stresses,double [] Ds)
        {
            var AA = 0.0;
            var BB = 0.0;
            var CC = 0.0;
            var dummy = 0.0;
            var XLamda = 0.0;
            AA = Math.Abs(Math.Pow(Ds[0], 2));
            BB = Ds[0] * (Stresses[0] - QH[13]);
            for (int i = 1; i < 6; i++)
            {
                BB = BB + (Ds[i] * (Stresses[i] - QH[13 + i])) / Math.Pow(PAR[3 + i], 2);
                AA = AA + Math.Pow((Ds[i] / PAR[3 + i]), 2);
            }
            if (AA == 0)
            {
                XLamda = 1;
            }
            CC = FPYE(PAR, QH, Stresses);
            dummy = Math.Pow(BB, 2) - AA * CC;
            XLamda = (-BB + Math.Sqrt(dummy)) / (AA);
            if (CC == 0)
            {
                XLamda = 1;
            }
            if (XLamda > 1)
            {
                XLamda = 1;
            }
            if (XLamda < 0)
            {
                XLamda = 0;
            }
            return XLamda;
        }
        #endregion
        #region CompSubroutines
        public void compsl (double [] PAR,double [] QH,double [] Stresses)
        {
            QH[13] = (1 - PAR[15]) * Stresses[0] + PAR[15] * QH[1];
            for (int i=1; i<6; i++)
            {
                QH[13 + i] = (1 - PAR[15]) * Stresses[i];
            }
            this.QH = QH;
        }
        public double compldot (double [] PAR,double[] de,double [] QH,double [] Stresses)
        {
            var Xn = 0.0;
            var Xd = 0.0;
            var H = 0.0;
            var Xldot = 0.0;
            Qgrad = compp(PAR, QH, Stresses);
            Xn = Quad(PAR, de, QH, Stresses, Qgrad);
            Xd = Quad(PAR, Qgrad, QH, Stresses, Qgrad);
            H = comph(PAR,QH,Stresses);
            Xd = Xd + H;
            if (Math.Abs(Xn) < Math.Pow(10, -15))
            {
                Xldot = 0.0;
            }
            else
            {
                if (Xd == 0)
                {
                    Xd = Math.Pow(10, -15);
                }
                Xldot = Xn / Xd;
            }
            return Xldot;
        }
        public double compyldot(double [] PAR,double [] QH,double [] Stresses,double [] ds,double [] de,double DA)
        {
            var yldot = 0.0;
            var Xn = 0.0;
            var Xd = 0.0;
            var Ratio = 0.0;
            Ratio = DA / QH[1];
            for (int i = 1; i < 6; i++)
                {
                  Xn = Xn + (Stresses[i]-QH[13+i]) * (ds[i]-Ratio*Stresses[i]) / (Math.Pow(PAR[3+i],2));
                  Xd = Xd + (Stresses[i] - QH[13 + i]) * (Stresses[i]) / (Math.Pow(PAR[3 + i], 2));
                }
            Xn=Xn+(Stresses[0]-QH[13]) * (ds[0]-Ratio* Stresses[0]);
            Xd=Xd+(Stresses[0] - QH[13]) * (Stresses[0]-QH[1]);
            Xd = PAR[15] * Math.Pow(QH[1], 2) - Xd;
            if (Xd == 0)
            {
                Xd = Math.Pow(10, -15);
            }
            yldot = Xn / Xd;
            if (Math.Abs(Xn) < Math.Pow(10, -15))
            {
                yldot = 0.0;
            }
            return yldot;
        }
        public double [] compp(double[] PAR, double[] QH, double[] Stresses)
        {
            Qgrad[0] = 2*(Stresses[0] - QH[13]);
            for (int i = 1; i < 6; i++)
            {
                Qgrad[i] = 2 * (Stresses[i] - QH[13 + i]) / Math.Pow(PAR[3+i],2);
            }
            return Qgrad;
        }
        public double[] comppbse(double[] PAR, double[] QH, double[] Stresses)
        {
            Pgrad[0] = 2 * (Stresses[0] - QH[1])*PAR[15];
            for (int i = 1; i < 6; i++)
            {
                Pgrad[i] = 2 *PAR[15]* (Stresses[i]) / Math.Pow(PAR[3 + i], 2);
            }
            return Pgrad;
        }
        public double compdotpdot (double [] Pgrad)
        {
            var Dpp = 0.0;
            for (int i=1; i<6; i++)
            {
                Dpp = Dpp + Pgrad[i]*Pgrad[i];
            }
            Dpp = Dpp / 1.5;
            Dpp = Math.Sqrt(Dpp);
            return Dpp;
        }
        public double [] conjugate(double[] PAR,double[] QH,double [] Stresses)
        {
            var SM = new double[6];
            SM[0] =QH[1]+(Stresses[0]-QH[13])/PAR[15];
            for (int i=1;i<6;i++)
            {
                SM[i] =(Stresses[i]-QH[13+i])/PAR[15];
            }
            return SM;
        }
        public double comphfact (double [] PAR,double [] QH,double [] Stresses,double [] SM)
        {
            var AA = 0.0;
            var BB = 0.0;
            var Halfterm1 = 0.0;
            var Halfterm2 = 0.0;
            var Hfact = 0.0;
            var Xd = 0.0;
            for (int i = 1; i<6;i++)
            {
                Halfterm1 =(SM[i]-Stresses[i])/PAR[3+i];
                Halfterm2 = (Stresses[i] - QH[7+i]) / PAR[3 + i];
                AA = AA + Halfterm1 * Halfterm1;
                BB = BB + Halfterm2 * Halfterm2;
            }
            AA = AA + Math.Pow((SM[0]-Stresses[0]), 2);
            BB = BB + Math.Pow((Stresses[0]-QH[7]), 2);
            if (Math.Abs(AA)<Math.Pow(10,-15))
            {
                Hfact = 0.0;
            }
            else
            {
                Qgrad = compp(PAR, QH, Stresses);
                Xd = Quad(PAR, Qgrad, QH, Stresses, Qgrad);
                if (Math.Abs(BB)==0)
                {
                    BB = Math.Pow(10, -15);
                }
                Hfact = PAR[16] * Xd * Math.Pow((AA / BB), PAR[17]);
            }
            return Hfact;
        }
        public double comph (double [] PAR,double [] QH,double [] Stresses)
        {
            var H = 0.0;
            var Hf1 = 0.0;
            var Hf2 = 0.0;
            var Dpart1 = 0.0;
            var Dpart2 = 0.0;
            var Dpart3 = 0.0;
            var Dexpart = 0.0;
            var Hpart1 = 0.0;
            var Hpart2 = 0.0;
            var Hpart3 = 0.0;
            var Dpp = 0.0;
            var SM = new double[6];
            var Hfact = 0.0;
            if (QH[3] == 1)
            {
                Qgrad = compp(PAR, QH, Stresses);
                Dpp = compdotpdot(Qgrad);
                Dexpart = Math.Exp(-Math.Abs(PAR[11] * QH[4] + PAR[12] * QH[5]));
                Dpart1 = Math.Abs(PAR[9]-PAR[10]+PAR[13] * QH[4] + PAR[14] * QH[5]);
                Dpart2 = Math.Abs(PAR[11] * Math.Abs(Qgrad[0]) + PAR[12] * Dpp);
                Dpart3 = Math.Abs(PAR[13] * Math.Abs(Qgrad[0]) + PAR[14] * Dpp);
                Hpart1 =(QH[2])*(Qgrad[0])*(Dpart1*Dexpart+PAR[10])/(PAR[2]-PAR[3]);
                Hpart2 = Math.Abs(Dpart3*Dexpart);
                Hpart3 =-Math.Abs(Dpart1*Dpart2*Dexpart);
                H =2*Stresses[0]*QH[6]*PAR[15]*(Hpart1+Hpart2+Hpart3);
                return H;
            }
            else
            {
                SM = conjugate(PAR, QH, Stresses);
                Hfact = comphfact(PAR, QH, Stresses, SM);
                Pgrad = comppbse(PAR, QH, SM);
                Dpp = compdotpdot(Pgrad);
                Dexpart = Math.Exp(-Math.Abs(PAR[11] * QH[4] + PAR[12] * QH[5]));
                Dpart1 = Math.Abs(PAR[9] - PAR[10] + PAR[13] * QH[4] + PAR[14] * QH[5]);
                Dpart2 = Math.Abs(PAR[11] * Math.Abs(Pgrad[0]) + PAR[12] * Dpp);
                Dpart3 = Math.Abs(PAR[13] * Math.Abs(Pgrad[0]) + PAR[14] * Dpp);
                Hpart1 = (QH[2]) * (Pgrad[0]) * (Dpart1 * Dexpart + PAR[10]) / (PAR[2] - PAR[3]);
                Hpart2 = Math.Abs(Dpart3 * Dexpart);
                Hpart3 = -Math.Abs(Dpart1 * Dpart2 * Dexpart);
                Hf1 =2*SM[0]* QH[6] * PAR[15] * (Hpart1 + Hpart2 + Hpart3);
                Hf2 =Hfact;
                QH[19] = Hf1;
                QH[20] = Hf2;
                H = Hf1 + Hf2;
                this.QH = QH;
                return H;
            }
        }
        #endregion
        #region Cstif/Compds/QuadSubroutines
        public IMatrix2D<double> cstif (double [] PAR,double [] QH,double [] Stresses)
        {
            var Xk = 0.0;
            var shear = 0.0;
            var dd = new double[6];
            var ee = new double[6];
            var omega = 0.0;
            var H = 0.0;
            Xk = QH[2] * Stresses[0] / PAR[3];
            shear = Xk * PAR[1];
            ConstitutiveMatrix[0, 0] = Xk;
            for (int i = 1; i<6;i++)
            {
                ConstitutiveMatrix[i, i] = shear;
            }
            if (QH[3] ==0)
            {
                return ConstitutiveMatrix;
            }
            else
            {
                Qgrad=compp(PAR, QH, Stresses);
                dd = compds(PAR, Qgrad, QH, Stresses);
                for (int i = 0; i< 6;i++)
                {
                    omega = omega +Qgrad[i]*dd[i];
                }
                H = comph(PAR, QH, Stresses);
                omega = omega + H;
                ee = compds(PAR, Qgrad, QH, Stresses);
                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 6; j++)
                    {
                        ConstitutiveMatrix[i, j] = ConstitutiveMatrix[i, j] - dd[i] * ee[j] / omega;
                    }
                return ConstitutiveMatrix;
            }
        }
        public double [] compds(double [] PAR,double [] de,double [] QH,double [] Stresses)
        {
            var ds = new double[6];
            var Xk = 0.0;
            var shear = 0.0;
            Xk = QH[2] * Stresses[0] / PAR[3];
            shear = Xk * PAR[1];
            ds[0] = Xk * de[0];
            for (int i = 1; i < 6; i++)
            {
                ds[i] = shear * de[i];
            }
            return ds;
        }
        public double Quad (double [] PAR,double [] de,double [] QH,double [] Stresses, double [] Qgrad)
        {
            var value = 0.0;
            var dd = new double[6];
            dd = compds(PAR, Qgrad, QH, Stresses);
            for (int i=0;i<6;i++)
            {
                value = value + Qgrad[i] * dd[i];
            }
            return value;
        }
        #endregion
    }
}
//George compile control and logic check