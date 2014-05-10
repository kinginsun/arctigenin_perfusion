;; 1. Based on: run2
;; 2. Description: 1-comp iv, linear elim
;; x1. Author: user
;; 3. Label:

$PROBLEM Fish PK analysis

$INPUT ID TIME DV AMT RATE MDV CMT EVID

$DATA fish1.csv IGNORE=@

$SUBROUTINES ADVAN6 TRANS1 TOL=3

$MODEL
COMP = (CENTRL)

$PK
Kout = THETA(1) * EXP(ETA(1))
Ke   = THETA(2) * EXP(ETA(2))
V    = THETA(3) * EXP(ETA(3))
S1   = V

$DES
IF(TIME<=30)THEN
DADT(1) = -Kout-Ke*A(1)
ELSE
DADT(1) = -Ke*A(1)
ENDIF

$ERROR
IPRED = F
Y = IPRED + EPS(1)

$THETA
(0, 14800) ; Kout
(0, 0.0471,1) ; Ke
(0, 1.51) ; V

$OMEGA
 0 FIX  ; IIV Kout
 0 FIX  ; IIV Ke
 0 FIX  ; IIV V

$SIGMA
 46500 ; Residual error

$EST METHOD=0 MAXEVAL=9999 NOABORT SIG=3 PRINT=1 POSTHOC
$COV

; Xpose
$TABLE ID TIME DV MDV CMT EVID IPRED ONEHEADER NOPRINT FILE=sdtab003.tab
$TABLE ID Kout Ke V ONEHEADER NOPRINT FIRSTONLY FILE=patab003.tab

