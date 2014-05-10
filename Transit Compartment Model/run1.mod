;; 1. Based on: 
;; 2. Description: 1-comp iv, linear elim
;; 3. Label:
;; x1. Author: user

$PROBLEM PK

$INPUT ID TIME DV AMT CMT MDV EVID

$DATA nm_001.csv ; IGNORE=@

$SUBROUTINES ADVAN1 TRANS2

$PK
CL = THETA(1) * EXP(ETA(1))
V  = THETA(2) * EXP(ETA(2))
S1 = V

$ERROR
IPRED = F
    W = SQRT(THETA(3)**2*IPRED**2 + THETA(4)**2)
    Y = IPRED + W*EPS(1)
 IRES = DV-IPRED
IWRES = IRES/W

$THETA
(0, 1)  ; CL
(0, 1)  ; V
(0, .1) ; Prop.RE (sd)
(0, 1)  ; Add.RE (sd)

$OMEGA
(0.1) ; IIV CL
(0.1) ; IIV V1

$SIGMA
1 FIX ; Residual error

$EST METHOD=1 INTER MAXEVAL=2000 NOABORT SIG=3 PRINT=1 POSTHOC
$COV

; Xpose
$TABLE ID TIME DV MDV EVID IPRED IWRES ONEHEADER NOPRINT FILE=sdtab001
$TABLE CL V ONEHEADER NOPRINT FIRSTONLY FILE=patab001
