for((i=15;i<45;i++))
do
mkdir rho-$i

for((j=140;j<180;j=j+2))
do
mkdir rho-$i/Temp-$j


echo "../../rho-0.$i 0 0 0 0 y y y n n ../../schedule_equil-$j y n dimer-react-$i-$j 1000 n n y dimer-txt-$i-$j n 1000 n y dimer-par-$i-$j n 1000 y y y y y y y 1000000" > start-$i-$j

echo "#$ -N e$i-$j" >equil$i-$j.sge
echo "#$ -P phaseamp" >>equil$i-$j.sge
echo "#$ -S /bin/sh" >>equil$i-$j.sge
echo "#$ -l h_rt=48:00:00" >>equil$i-$j.sge
echo "/bin/echo Running on host: `hostname`.">>equil$i-$j.sge
echo "/bin/echo Starting on: `date`" >>equil$i-$j.sge
echo "cd rho-$i/Temp-$j" >>equil$i-$j.sge
echo "/bin/echo In directory: `pwd`" >>equil$i-$j.sge
echo "../../dmd-react/dmd.out < ../../start-$i-$j> ./out" >> equil$i-$j.sge
echo "cd ../../." >> equil$i-$j.sge
qsub equil$i-$j.sge

done; 

done



for((i=80;i<100;i++))
do
mkdir rho-$i

for((j=100;j<140;j=j+2))
do
mkdir rho-$i/Temp-$j


echo "../../rho-0.$i 0 0 0 0 y y y n n ../../schedule_equil-$j y n dimer-react-$i-$j 1000 n n y dimer-txt-$i-$j n 1000 n y dimer-par-$i-$j n 1000 y y y y y y y 1000000" > start-$i-$j

echo "#$ -N e$i-$j" >equil$i-$j.sge
echo "#$ -P phaseamp" >>equil$i-$j.sge
echo "#$ -S /bin/sh" >>equil$i-$j.sge
echo "#$ -l h_rt=48:00:00" >>equil$i-$j.sge
echo "/bin/echo Running on host: `hostname`.">>equil$i-$j.sge
echo "/bin/echo Starting on: `date`" >>equil$i-$j.sge
echo "cd rho-$i/Temp-$j" >>equil$i-$j.sge
echo "/bin/echo In directory: `pwd`" >>equil$i-$j.sge
echo "../../dmd-react/dmd.out < ../../start-$i-$j> ./out" >> equil$i-$j.sge
echo "cd ../../." >> equil$i-$j.sge
qsub equil$i-$j.sge

done; 

done


