#!bin/zsh

NomalSolver()
{
# $1 -> Matrix $2 -> solver
  rm -f ./output/*
  ./Solver -m $1 -S $2  >> ../$DIR/$MAT/$1-$2.log
  cp ./output/${2^^}_his.txt ../$DIR/$MAT/$1-$2.txt
  echo "-------- Matrix [$1] solver [$2] Done ------"
  gnuplot -e "INPUTNAME='../$DIR/$MAT/$1-$2.txt'; OUTPUTNAME='$1-$2.eps'; TITLE='$2'" ./plot.gnp
  mv ./$1-$2.eps ../$DIR/$MAT/EPS/
  echo "-------- Plot Done ------"
}

VPSolver()
{
  for i in `seq 1 $4`
  do
# $1 -> Matrix $2 -> outer sovler $3 -> inner solver $4 max innerloop 
    rm -f ./output/*
    ./Solver -m $1 -S $2 -s $3 -k $i -l 10 -R 1000 >> ../$DIR/$MAT/$1-$2-$3-k$i.log
    cp ./output/${2^^}_his.txt ../$DIR/$MAT/$1-$2-$3-k$i.txt
    echo "-------- Matrix [$1] OuterSolver [$2] InnerSolver [$3] in k [$i] Done ------"
  done
  gnuplot -e "OUTPUTNAME='../$DIR/$MAT/EPS/$1-$2-$3-k.eps'; P1='../$DIR/$MAT/$1-$2-$3-k1.txt'; N1='$2-$3-k1'; P2='../$DIR/$MAT/$1-$2-$3-k2.txt'; N2='$2-$3-k2'; P3='../$DIR/$MAT/$1-$2-$3-k3.txt'; N3='$2-$3-k3'; P4='../$DIR/$MAT/$1-$2-$3-k4.txt'; N4='$2-$3-k4'" ./plot_mult.gnp
  echo "-------- Plot Done ------"
}

OneLine()
{
  echo "-------- Plot Oneline ------"
  gnuplot -e "OUTPUTNAME='../$DIR/$MAT/EPS/$1.eps'; P1='../$DIR/$MAT/$1-$2.txt'; N1='$2'; P2='../$DIR/$MAT/$1-$3.txt'; N2='$3'; P3='../$DIR/$MAT/$1-$4.txt'; N3='$4'; P4='../$DIR/$MAT/$1-$5.txt'; N4='$5'; P5='../$DIR/$MAT/$1-$6.txt'; N5='$6'" ./plot_onetime.gnp
  echo "-------- Plot Done ------"
}

echo "===================start========================"
if [ $# -ne 1 ]; then
  echo "指定された引数は$#個です。"
  echo "実行するには1個の引数が必要です。"
  exit 1
fi

DIR="16-9-20-L10"
MAT=$1

DIRFULL=../$DIR/$MAT/EPS/
if [ ! -e $DIRFULL ]; then
  mkdir -p $DIRFULL
  echo "mkdir"
fi

#NomalSolver $MAT "cg"
#NomalSolver $MAT "cr"
#NomalSolver $MAT "gcr"
#NomalSolver $MAT "gmres"
#NomalSolver $MAT "kskipcg"
#

#VPSolver $MAT "vpcg" "kskipcg" "4"
#VPSolver $MAT "vpgcr" "kskipcg" "4"
VPSolver $MAT "vpgmres" "kskipcg" "4"

#OneLine $MAT "cg" "cr" "gcr" "gmres" "kskipcg" 

echo "===================over========================="

