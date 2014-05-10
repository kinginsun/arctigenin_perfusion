;; 1. Based on: run106
;; 2. Description: Test Adv6 CMT2
;; x1. Author: user

$PROBLEM  Test 2 Single Comp- ADVAN6
$DATA   AVD6test2.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV
$SUBROUTINE ADVAN6 TOL=6
$MODEL
COMP=(CENTRAL,DEFOBS,DEFDOSE)
COMP=(PP)
$PK
K10= THETA(1)*EXP(ETA(1))
V1= THETA(2)*EXP(ETA(2))
S1= V1

$DES
DADT(1)=-K10*A(1)

$ERROR
IPRED = A(1)/V1
W = 1;SQRT(THETA(3)**2*IPRED**2 + THETA(4)**2)
Y = IPRED + W*EPS(1)
IRES = DV-IPRED
IWRES = IRES/W

$THETA
(0, 1)  ; K10
(0, 5)  ; V1

$OMEGA
(0.01) ; IIV K10
(0.01) ; IIV V1


$SIGMA
0.01 ; Residual error

;；$SIM (12345) (54321) ONLYSIM
$EST METHOD=0 MAXEVAL=9999 NOABORT SIG=6 PRINT=5 POSTHOC
$COV
; ; ; Xpose
$TABLE ID TIME DV MDV IPRED IWRES ONEHEADER NOPRINT FILE=sdtab107.tab
$TABLE ID K10 V1 ETA1 ETA2 ONEHEADER NOPRINT FILE=patab107.tab
