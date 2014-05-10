;; 1. Based on: run1
;; 2. Description: 1-comp iv, linear elim
;; x1. Author: user
;; 3. Label:

$PROBLEM Fish PK analysis

$INPUT ID TIME DV AMT RATE MDV CMT EVID

$DATA fish.csv IGNORE=@

$SUBROUTINES ADVAN6 TRANS1 TOL=3

$MODEL
COMP = (ABSORP)
COMP = (CENTRL)

$PK
Kout = THETA(1) * EXP(ETA(1))
Ka   = THETA(2) * EXP(ETA(2))
Ke   = THETA(3) * EXP(ETA(3))
V    = THETA(4) * EXP(ETA(4))
S2   = V

$DES
DADT(1) = -Kout-Ka*A(1)
DADT(2) = Ka*A(1)-Ke*A(2)

$ERROR
IPRED = F
Y = IPRED + EPS(1)

$THETA
(0, 1400)  	; Kout
(0, 5)  	; Ka
(0, .05, 1) 	; Ke
(0, 100)  	; V

$OMEGA
0 FIX ; IIV Kout
0 FIX ; IIV Ka
0 FIX ; IIV Ke
0 FIX ; IIV V

$SIGMA
10 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=3 PRINT=1 POSTHOC
$COV

; Xpose
$TABLE ID TIME DV MDV CMT EVID IPRED ONEHEADER NOPRINT FILE=sdtab002.tab
$TABLE ID Kout Ka Ke V ONEHEADER NOPRINT FIRSTONLY FILE=patab002.tab

