BASEDIR=` pwd `

echo running $1
qsub -v infile=`pwd`/$1,addinfiles=`pwd`/pndriftHCT.m:`pwd`/EquilibratePNHCT.m:`pwd`/hstruct.m:`pwd`/p3hTPCBM.csv:`pwd`/readlayersCT.m:`pwd`/DriftanalyseCT.m:`pwd`/v2struct.m:`pwd`/pnParamsHCT.m:`pwd`/symmetricize.m:`pwd`/tpvfit.m:`pwd`/importdata1.m  Matlab.sh
