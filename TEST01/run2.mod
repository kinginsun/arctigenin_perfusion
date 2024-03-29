;; 1. Based on: run1
;; 2. Description: 1-comp iv, linear elim
;; x1. Author: user
;; 3. Label:

$PROBLEM PK

$INPUT      ID TIME CP=DV DOSE=AMT MDV
$DATA       test_data.txt 

$SUBROUTINES  ADVAN2

$PK
;THETA(1)=MEAN ELIMINATION RATE CONSTANT K (1/HR)
;THETA(2)=VOLUME OF DISTRIBUTION
;THETA(3)=MEAN ABSORPTION RATE CONSTANT KA (1/HR)
CALLFL=1
K=THETA(1)+ETA(1) 
KA=THETA(2)+ETA(2) 
S2= THETA(3)+ETA(3)
;CL=K*V ;


$THETA
(0.001, 0.625,5) FIX ;
(0, 1.61,0) FIX ;
(0, 0.458,0) FIX ;

$OMEGA BLOCK(3) 0.276  ;VARIANCE OF ETA

$ERROR
Y=F+EPS(1)

$SIGMA 0.909 FIX  ;VARIANCE OF EPS  

$SIM (12345) (54321) ONLYSIM
; $EST METHOD=1 INTER MAXEVAL=2000 NOABORT SIG=3 PRINT=1 POSTHOC
; $COV
; ; ; Xpose
$TABLE ID TIME DV MDV PRED WRES ONEHEADER NOPRINT FILE=sdtab002
$TABLE KA K ONEHEADER NOPRINT FIRSTONLY FILE=patab002

