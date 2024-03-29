;; 1. Based on: run109
;; 2. Description: 1-comp iv, MM elim
;; x1. Author: user
;; 3. Label:
$PROBLEM Simulation of population data and weight
$INPUT ID TIME WT AMT DV
$DATA SIMORIG.txt IGNORE=#
$SUBROUTINE ADVAN1
$PK
SIMWT=70+70*ETA(3)
CL=THETA(1)*SIMWT*EXP(ETA(1))
V=THETA(2)*EXP(ETA(2))
K=CL/V
S1=V
$ERROR
Y=F+F*EPS(1)
$THETA  .0625  10
$OMEGA .09 .05 .04
$SIGMA .01
$SIMULATION (9215690) ONLYSIM SUBPROBLEMS=10
$TABLE ID TIME SIMWT AMT NOPRINT FILE=SIMDATA2 NOHEADER

