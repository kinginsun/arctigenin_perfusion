;; 1. Based on: 
;; 2. Description: PK-PD effect compartment
;; 3. Label:
;; x1. Author: user

$PROBLEM PK-PD

$INPUT ID TIME DV AMT CMT MDV EVID

$DATA nm_001.csv ; IGNORE=@

$SUBROUTINES ADVAN6 TRANS1 TOL=5

$MODEL
  COMP=(PK)
  COMP=(EFF)

$PK
  CL   = THETA(1) * EXP(ETA(1))
  V1   = THETA(2) * EXP(ETA(2))
  KE   = CL/V1
  S1   = V1
  KE0  = THETA(3) * EXP(ETA(3))
  EMAX = THETA(4)
  EC50 = THETA(5) * EXP(ETA(4))
  E0   = THETA(6) * EXP(ETA(5))

$DES
  DADT(1) = -KE * A(1)
  DADT(2) = KE0*(A(1)- A(2))

$ERROR
IF (CMT.EQ.1) THEN
  IPRED = F
    W   = SQRT(THETA(7)**2*IPRED**2 + THETA(8)**2)
    Y   = IPRED + W*EPS(1)
ELSE
  IPRED = E0 + (EMAX*A(2))/(EC50+A(2))
  W     = THETA(9)
  Y     = IPRED + W*EPS(2)
ENDIF
W = 1
IRES = DV-IPRED
IWRES = IRES/W

$THETA
(0,1) ; CL
(0,1) ; V
(0,1) ; KE0
(1)   ; EMAX
(1)   ; EC50
(1)   ; E0
(0, .1) ; Prop.RE (sd)
(0, 1)  ; Add.RE (sd)
(0, 1)  ; Add.RE PD (sd)

$OMEGA
(0.1) ; IIV CL
(0.1) ; IIV V
(0.1) ; IIV KE0
(0.1) ; IIV EC50
(0.1) ; IIV E0

$SIGMA
1 FIX
1 FIX ; Additive error

$EST METHOD=1 INTER MAXEVAL=2000 NOABORT SIG=3 PRINT=1 POSTHOC
$COV

$TABLE ID TIME DV MDV EVID IPRED IWRES ONEHEADER NOPRINT FILE=sdtab008
$TABLE CL V1 E0 EMAX EC50 KE0 FIRSTONLY ONEHEADER NOPRINT FILE=patab007
