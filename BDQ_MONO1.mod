$PROBLEM 19-state PD model - PKPD model
;; 1. Based on: BDQ_MONO
;; 2. Description: bdq model with changes, 8 cmt
;; x1. Author: peter
 
 
$INPUT ID	DV	MDV	EVID	AMT	TIME	CMT	DINT	DOSE	DRUG	IBMAX	RATE	MONO	ICL	IV2	IF	IF2	HIV	CD4	CAVITY	WT	FEMALE	IV11	IKA	IMTT	ICLL	IF11
 
$DATA EXTRAP_B300_2.csv IGNORE=@ IGNORE=(DRUG.NE.1) IGNORE=(ID.EQ.9025) IGNORE=(ID.EQ.9063) IGNORE=(ID.EQ.8023) IGNORE=(ID.EQ.9025) IGNORE=(ID.EQ.9063) IGNORE=(ID.EQ.8023) IGNORE=(ID.EQ.1005) IGNORE=(ID.EQ.1018) IGNORE=(ID.EQ.2008) IGNORE=(ID.EQ.2001) IGNORE=(ID.EQ.8025) IGNORE=(ID.EQ.8005) IGNORE=(ID.EQ.8023) IGNORE=(ID.EQ.8042) IGNORE=(ID.EQ.8055) IGNORE=(ID.EQ.8062) IGNORE=(ID.EQ.9025) IGNORE=(ID.EQ.9037) IGNORE=(ID.EQ.9049) IGNORE=(ID.EQ.9068)  
IGNORE=(ID.EQ.9071) 
 
$SUBROUTINE ADVAN13 TOL=9
 
$MODEL NCOMP=8
 
$OMEGA
0 FIX 
0 FIX 
0 FIX 
0 FIX 
 
$PK
 
IF (NEWIND.NE.2) THEN
 
LAUC = 0
 
LTM = TIME
 
LAVE = 0
 
ENDIF
 
rF = THETA(1) ; knet growth death rate - fast state (/h)
 
rS = THETA(2) ; knet growth death rate - slow state (/h)
 
TVS0 = THETA(3)
 
KE0 = THETA(19)*EXP(ETA(1))
IF (DOSE.EQ.200) KE0 = THETA(23)*EXP(ETA(2))
IF (DOSE.EQ.300) KE0 = THETA(24)*EXP(ETA(3))
IF (DOSE.EQ.400) KE0 = THETA(25)*EXP(ETA(4))
 
KE1 = THETA(20)
 
BMAX = (10**IBMAX)
 
A_0(1) = (10**TVS0)
 
TVCL = THETA(4)

CL = ICL
 
TVV4 = THETA(5)
 
V4=IV2
 
Q3 = THETA(6)
V5 = THETA(7)
Q4 = THETA(9)
V6 = THETA(10)
Q5 = THETA(11)
V7 = THETA(12)
F3=IF
F4=IF2
D3 = THETA(13)
ALAG4 = D3 + THETA(14)
D4 = THETA(15)
ALAG3 = THETA(8)
K34=1000
K20 = CL/V4
K45 = Q3/V4
 
K54 = Q3/V5
 
K46 = Q4/V4
 
K64 = Q4/V6
 
K47 = Q5/V4
 
K74 = Q5/V7
 
S4 = V4
 
EMAX_F = THETA(16)
 
UNBOUND_MICE = 0.01
UNBOUND_HUMAN = 0.01
 
UNBOUND_MICE_PA = 0.08
UNBOUND_HUMAN_PA = 0.14
 
; MONOTHERAPY POTENCY 1005/1008
 
IC50_F = THETA(17)
 
IC50_S = THETA(18)
 
; DDI EFFECT
 
PA_EFF_F = 1
 
PA_EFF_S = 1
 
PZ_EFF_F = 1
 
PZ_EFF_S = 1
 
IF (DRUG.EQ.12) PA_EFF_F = THETA(21)
IF (DRUG.EQ.12) PA_EFF_S = THETA(21)
 
NEW_IC50_F = (IC50_F*(UNBOUND_MICE/UNBOUND_HUMAN)*PZ_EFF_F*PA_EFF_F)
NEW_IC50_S = (IC50_S*(UNBOUND_MICE/UNBOUND_HUMAN)*PZ_EFF_S*PA_EFF_S)
NEW_EMAX = (EMAX_F) ; ASSUME EMAX AGAINST F IS SAME AS S POPULATION
 
AUC = A(8)
 
DTIM = TIME - LTM ; Time difference
 
IF (DTIM.EQ.0.AND.NEWIND.EQ.2) AUC = LAUC
 
IF (DTIM.EQ.0.AND.NEWIND.EQ.2) LAUC = 0
 
DAUC = AUC - LAUC
 
IF (DTIM.EQ.0) THEN
 
AVECONC = LAVE
 
ELSE
 
AVECONC = DAUC/DINT
 
ENDIF
 
LAUC = AUC
 
LTM = TIME
 
LAVE = AVECONC
 
RATIO=THETA(22)
 
EFF_FBUGS = (NEW_EMAX*AVECONC/(NEW_IC50_F+AVECONC))
 
EFF_SBUGS = (NEW_EMAX*AVECONC/(NEW_IC50_S+AVECONC))
 
$DES
 
kFS = (rF*(A(1)+A(2)))/BMAX
 
kSF = kFS/RATIO
 
GROWTHFUNC_F = rF*(1-(A(1)+A(2))/BMAX)
 
GROWTHFUNC_S = rS*(1-(A(1)+A(2))/BMAX)
 
IF (GROWTHFUNC_F.LT.0) GROWTHFUNC_F = 0
 
IF (GROWTHFUNC_S.LT.0) GROWTHFUNC_S = 0
 
GROWTH_F = (GROWTHFUNC_F-EFF_FBUGS)
GROWTH_S = (GROWTHFUNC_S-EFF_SBUGS)
 
DADT(1) = GROWTH_F*A(1)-kFS*A(1)+kSF*A(2)
DADT(2) = GROWTH_S*A(2)+kFS*A(1)-kSF*A(2)
DADT(3) = -K34*A(3)
DADT(4) = K34*A(3)-K20*A(4)- K45*A(4)-K46*A(4)-K47* A(4)+K54* A(5)+ K64* A(6)+K74*A(7)
DADT(5) = K45* A(4)- K54* A(5)
DADT(6) = K46* A(4)- K64* A(6)
DADT(7) = K47* A(4)- K74* A(7)

 
;; Cumulative AUC
 
DADT(8) = A(4)/V4
 
$ERROR
 
FBUGS = 0
 
IF(A(1).GE.1) FBUGS = A(1)
 
SBUGS = 0
 
IF(A(2).GE.1) SBUGS = A(2)
 
;IPRED = 0
 
;IF (FBUGS.GE.1.OR.SBUGS.GE.1) IPRED = LOG10(FBUGS+SBUGS)

IPRED = F

Y= F + ERR(1)
 
;;add BDQ thetas
 
$THETA
(0, 0.0272) FIX ; 1 rF - from BALB/c
(0, 0.00068) FIX ; 2 rS - from BALB/c
(0, 2) FIX ; 3 inoculum in log10 scale - fixed
(0.01, 2.77) FIX ; POPCL (4)
(0.01, 163) FIX ; POPV2 (5)
(0.01, 11.7) FIX ; POPQ3 (6)
(0.01, 178) FIX ; POPV3 (7)
(0.01, 0.916) FIX ; ALAG (8)
(0.01, 8.03) FIX ; POPQ4 (9)
(0.01, 2990) FIX ; POPV4 (10)
(0.01, 3.41) FIX ; POPQ5 (11)
(0.01, 7360) FIX ; POPV5 (12)
(0.01, 2.21) FIX ; D1B (13)
(0.01, 1.47) FIX ; TLAGB (14)
(0.01, 1.47) FIX ; D2B (15)
(0.0671) FIX ; 16 EMAX of RIF
(0.192) FIX ; 17 IC50 F for BDQ (mg/L)
(3.04) FIX ; 18 IC50 S for BDQ (mg/L)
(0, 0.0181) ; 19 KE0
(0, 0.0072) FIX ; 20 KE1
(1.72) FIX ; 21 BPA
(0.7) FIX ; 22 Ratio
(0, 0.0567) ; 23 KE0
(0, 0.0406) ; 24 KE0
(0, 0.0457) ; 25 KE0
 
$SIGMA
0.409
 
$TABLE ID TIME DV IPRED HIV CD4 MDV AVECONC FBUGS SBUGS CMT AUC DOSE DRUG NOPRINT ONEHEADER FILE=sdtabBDQ13

$SIM (20030521 NORMAL NEW) ONLYSIM SUBPROBLEMS=100
 
;$ESTIMATION MAXEVAL=2000 METHOD=1 INTER NOABORT POSTHOC NSIG=3 SIGL=9
;$COVARIANCE
