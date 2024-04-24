1;; 1. Based on: run727
;; 2. Description: Svensson pk model - approved dose
;; x1. Author: user
$PROBLEM    TMC207
$INPUT      ID WOMEN RACE DV TIME EVID CMT MDV AMT STUD WT ALB AGE OCC
            FLAG
$DATA      COSVEN4_WHO.csv IGNORE=@
$SUBROUTINE ADVAN13 TOL=6
$MODEL      NCOMPARTMENTS=8 COMP=(DEPOT DEFDOSE) COMP=(BDQC)
            COMP=(BDQPERI1) COMP=(BDQPERI2) COMP=(M2) COMP=(TRANSI1)
            COMP=(TRANSI2) COMP=(ALBUMIN)
$PK

;--- Model for Albumin

TVX0 = THETA(1)

TVXSS = THETA(2)

TVREP = THETA(3)

BSVX0 = ETA(1)

BSVXSS = ETA(2)

BSVREP = ETA(3)

SHPX0 = THETA(4)

PHIX0 = EXP(BSVX0)

PHI2X0 =(PHIX0**SHPX0-1)/SHPX0 ; Box-cox transformation

SHPXSS = THETA(5)

PHIXSS = EXP(BSVXSS)

PHI2XSS =(PHIXSS**SHPXSS-1)/SHPXSS ; Box-cox transformation

X0 = TVX0*EXP(PHI2X0)

XSS = TVXSS*EXP(PHI2XSS)

REP = TVREP*EXP(BSVREP)

HL = LOG(2)/REP

A_0(8) = X0

;--- Model for WT

TVWT0 = THETA(6)

TVWT120 = THETA(7)

BSVWT0 = ETA(4)

BSVWT120 = ETA(5)

SHPWT120 = THETA(8)

PHIWT120 = EXP(BSVWT120)

PHI2WT120 = (PHIWT120**SHPWT120-1)/SHPWT120 ; Box-cox transformation

WT0 = TVWT0*EXP(BSVWT0)

WT120 = TVWT120*EXP(PHI2WT120)

SLOPE = (WT120 - WT0)/(120*7*24) ; TIME in hours, 120 weeks

;--- BDQ and M2 PK

;--- Typical values fixed effects

TVMAT = THETA(9)

TVFR = THETA(10)

TVCL= THETA(11)

TVV= THETA(12)

TVQ1 = THETA(13)

TVVP1 = THETA(14)

TVQ2 = THETA(15)

TVVP2 = THETA(16)

TVCLM2 = THETA(17)

TVVM2 = THETA(18)

;--- Typical values variability

BOVF = ETA(6)

IF(OCC.EQ.2) BOVF = ETA(7)

BOVMAT = ETA(8)

IF(OCC.EQ.2) BOVMAT = ETA(9)

BSVF = ETA(10)

BSVCL= ETA(11)

BSVCLM2= ETA(12)

BSVV= ETA(13)

BSVQ1= ETA(14)

BSVVM2= ETA(15)


;--- Mechanistic

; Allometric scaling and albumin effects coded in $DES and $ERROR since they are time changing

;--- Empiric

; Effect of Black race on CL

BLACK = 0

IF(RACE.EQ.2) BLACK=1

BLACKCL = 1 + BLACK*THETA(23)

; AGE on CL

AGECL = 1 + (32-AGE)*(THETA(24))

;--- Parameters

PHI = LOG(TVMAT/(1-TVMAT))+BOVMAT

MAT = 6*EXP(PHI)/(EXP(PHI)+1)

; Mean absorption time, overall time for both delay and 90% complete absorption, logit transformed

; to retain constraines even with BOV in MAT

FR = TVFR ; Fraction of MAT that is delay

MTT = MAT*FR

KAHL= MAT*(1-FR)/3.3

KA = LOG(2)/KAHL

KTR= 2/MTT

F1= 1.8002*EXP(BOVF+BSVF)

; AMT in mg in input file, MW TMC207 555.5 g/mol, DV as nmol/mL = μmol/L -->

;(AMT/1000)/(555.5)*1000000 = AMT*1.8002 μmol

CLB = TVCL*BLACKCL*AGECL*EXP(BSVCL)

VB = TVV*EXP(BSVV)

Q1B= TVQ1*EXP(BSVQ1)

VP1B = TVVP1

Q2B = TVQ2

VP2B= TVVP2

CLM2B = TVCLM2*BLACKCL*AGECL*EXP(BSVCLM2)

VM2B= TVVM2*EXP(BSVVM2)

$DES

; --- Albumin

DADT(8) = LOG(2)/(HL*7*24)*A(8)*(1- A(8)/XSS) ; Time in hours, unit of half life: weeks

; --- Body weight

WTTIME = WT0 + T*SLOPE ; Predicted individual WT at time T

; --- Time varying covariates

ALBRELI = A(8)/XSS

; Time varying individual albumin relative individual value albumin at SS

COVALBI = (ALBRELI)**THETA(22) ; Albumin effect on hepatic function --> CL

FM = (ALBRELI)**(-THETA(22)) ; Albumin effect on hepatic function --> fm

ALLCL = (WTTIME/70)**THETA(19) ; Allometric scaling CL/Q

ALLV= (WTTIME/70)**THETA(20) ; Allometric scaling V/VP

; --- BDQ and M2 PK

CL = CLB*COVALBI*ALLCL

V = VB*ALLV

Q1 = Q1B*ALLCL

VP1 = VP1B*ALLV

Q2 = Q2B*ALLCL

VP2 = VP2B*ALLV

CLM2 = CLM2B/FM*COVALBI*ALLCL

VM2 = VM2B/FM*ALLV

DADT(1)= -KTR*A(1)

DADT(2) = A(7)*KA - A(2)*CL/V - A(2)*Q1/V + A(3)*Q1/VP1 - A(2)*Q2/V + A(4)*Q2/VP2 ; BDQ

DADT(3) = A(2)*Q1/V - A(3)*Q1/VP1

DADT(4) = A(2)*Q2/V - A(4)*Q2/VP2

DADT(5) = A(2)*CL/V - A(5)*CLM2/VM2 ; M2

DADT(6) = A(1)*KTR - A(6)*KTR ; transit1

DADT(7) = A(6)*KTR - A(7)*KA ; transit2

$ERROR

; --- Body weight

WTTIMEE = WT0 + TIME*SLOPE

; --- Time varying covariates

FURELIP = THETA(2)/A(8)

; Time varying individual fraction unbound relative typical population value at SS (inversely

;proportional to albumin)

COVFUIP = (FURELIP)**THETA(21) ; Fu effect on V

ALBRELIE = A(8)/XSS

; Time varying individual albumin relative individual value albumin at SS

FME = (ALBRELIE)**(-THETA(22)) ; Albumin effect on hepatic function --> fm

ALLCLE = (WTTIMEE/70)**THETA(19) ; Allometric scaling CL/Q

ALLVE = (WTTIMEE/70)**THETA(20) ; Allometric scaling V/VP

; --- BDQ and M2 PK

VE = VB*ALLVE*COVFUIP

VM2E = VM2B/FME*ALLVE*COVFUIP

; FLAG 1 = BDQ PK, 2= M2 PK, 3= Albumin, 4= Body weight

DEL= 1E-12

IPRED=LOG(A(2)/VE+DEL)

IF(FLAG.EQ.2) IPRED=LOG(A(5)/VM2E+DEL)

IF(FLAG.EQ.3) IPRED = A(8)

IF(FLAG.EQ.4) IPRED = WTTIMEE

; Error additive on log scale for PK, proportional for albumin and weight

BSVRUV1= ETA(16)

BSVRUV2= ETA(17)

W = 1*EXP(BSVRUV1)

IF(FLAG.EQ.2) W = 1*EXP(BSVRUV2)

IF(FLAG.EQ.3) W = IPRED

IF(FLAG.EQ.4) W = IPRED

IF(W.EQ.0) W=1

IRES=DV-IPRED

IWRES=IRES/W

Y = IPRED + W*EPS(3)

IF(FLAG.EQ.2) Y = IPRED + W*EPS(4)

IF(FLAG.EQ.3) Y = IPRED + W*EPS(1)

IF(FLAG.EQ.4) Y = IPRED + W*EPS(2)

$THETA  (0,3.64865) FIX ; 1 X0 g/dl
 (0,4.04068) FIX ; 2 Xss g/dl
 (0,0.03399109) FIX ; 3 Rate constant return to normal 1/week
 -2.43936 FIX ; 4 shape factor BSVX0 boxcox
 -5.37749 FIX ; 5 shape factor BSVXSS boxcox
 (0,56.6371) FIX ; 6 WT0 kg
 (0,62.6425) FIX ; 7 WT120 kg
 -0.416034 FIX ; 8 Boxcox BSV WT120
 (1E-06,0.662045,1) FIX ; 9 MAT
 (1E-06,0.466443,1) FIX ; 10 FR
 (1E-06,2.61592) FIX ; 11 CL
 (1E-06,198.34) FIX ; 12 V
 (1E-06,3.658) FIX ; 13 Q1
 (1E-06,8549.06) FIX ; 14 VP1
 (1E-06,7.33504) FIX ; 15 Q2
 (1E-06,2690.91) FIX ; 16 VP2
 (1E-06,10.0496) FIX ; 17 CLM2
 (1E-06,2203.71) FIX ; 18 VM2
 (1E-06,0.180912) FIX ; 19 Allometric scaling baseline CL
 (1E-06,1) FIX ; 20 Allometric scaling V
 (1E-06,1) FIX ; 21 Time varying FU on BDQ+M2 disp
 (-10,1.64021,10) ; 22 Individual time varying effect of ALB BDQ+M2 CL
 (0,0.838743) ; 23 Effect of black race on CL/CM2
 (0,0.00880756,10) ; 24 AGE effect on CL/CLM2
; --- Albumin
; --- Body weight
; --- BDQ and M2 PK
; --- Covariates
; --- Albumin and body weight
$OMEGA  BLOCK(5) FIX
 0.0253964  ;    1 BSVX0
 0.00787165 0.00965573  ;   2 BSVXSS
 -0.0831085 0.00887115 0.979314  ;   3 BSVREP
 0.0117627 0.00499615 -0.0379866 0.0421083  ;   4 BSVWT0
 0.00559706 0.00649501 -0.00226264 0.0369625 0.0494914  ; 5 BSVRWT120
; --- BDQ and M2 PK
$OMEGA  BLOCK(1) FIX
 0  ;    6 BOV F
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) FIX
 0  ;  8 BOV MAT
$OMEGA  BLOCK(1) SAME
$OMEGA  0.0803271  ;   10 BSV F
$OMEGA  BLOCK(2)
 0.152776  ;  11 BSV CL
 0.13486 0.212124  ; 12 BSV CLM2
$OMEGA  0.171909  ;   13 BSV V
 0.181182  ;  14 BSV Q1
 0.150223  ; 15 BSV VM2
$OMEGA  BLOCK(2)
 0.05392  ; 16 BSV RUVBDQ
 0.0295137 0.0522882  ; 17 BSV RUVM2
$SIGMA  0.00500974  ; 1 Prop error ALB
 0.00114573  ; 2 Prop error WT
$SIGMA  BLOCK(2)
 0.0518161  ; 3 Prop error TMC
 0.0189319 0.0366836  ; 4 Prop error M2
$SIMULATION (12345) (54321) ONLYSIM NSUBPROBLEMS=5
;

;$ESTIMATION MAXEVAL=2000 METHOD=1 INTER NOABORT NSIG=2

;$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 PRINT=1

;$COVARIANCE MATRIX=R UNCOND PRINT=E
$TABLE      ID TIME IPRED CWRES DV MDV CMT FILE=sdtab3

