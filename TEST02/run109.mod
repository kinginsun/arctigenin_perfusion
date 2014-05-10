;; 1. Based on: run108
;; 2. Description: 1-comp iv, MM elim
;; x1. Author: user
;; 3. Label:

$PROBLEM Simulation of population data
$INPUT ID TIME WT AMT DV
$DATA SIMORIG.txt IGNORE=#
$SUBROUTINE ADVAN1
$PK
CL=THETA(1)*EXP(ETA(1))
V=THETA(2)*EXP(ETA(2))
K=CL/V
S1=V
$ERROR
Y=F+F*EPS(1)
$THETA  .0625  10
$OMEGA .09 .05
$SIGMA .01
$SIMULATION (9215690) ; seed 1-7 digits
$TABLE ID TIME WT AMT NOPRINT FILE=SIMDATA1 NOHEADER


