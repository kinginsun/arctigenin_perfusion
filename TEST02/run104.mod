;; 1. Based on: run103
;; 2. Description: A SIMPLE PHARMACODYNAMIC MODEL SD=5
;; x1. Author: user

$PROBLEM A SIMPLE PHARMACODYNAMIC MODEL
$INPUT   ID TIME CP DV
$DATA    SinglePD.csv IGNORE=@
$PRED
EMAX=THETA(1)+ETA(1)
C50=THETA(2)+ETA(2)
E=EMAX*CP/(C50+CP)
Y=E+ERR(1)
$THETA   100   20
$OMEGA   0 FIX   0 FIX
$SIGMA   4
$ESTIMATION

