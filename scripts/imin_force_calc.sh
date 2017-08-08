MFP=GGGGG_vac.30
num=49-53

sander -O -i LJ2_q1.00.frc.nb.in -p $MFP.no_q.no1-4.prmtop -c $MFP.run50.rst7 -y $MFP.run$num.nc -x $MFP.run$num.vdw.nc -r $MFP.run$num.vdw.rst7 -inf $MFP.run$num.vdw.inf -o $MFP.run$num.vdw.log