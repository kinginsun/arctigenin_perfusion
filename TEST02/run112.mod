;; 1. Based on: run111
;; 2. Description: A SIMPLE PHARMACODYNAMIC MODEL
;; x1. Author: user

$PROBLEM A SIMPLE PHARMACODYNAMIC MODEL
$INPUT   ID CP DV
$DATA    SinglePD.csv IGNORE=@
$PRED
EMAX=THETA(1)
C50=THETA(2)
E=EMAX*CP/(C50+CP)
Y=E+ERR(1)
$THETA   100   20
$OMEGA   4
$ESTIMATION
