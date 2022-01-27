for((i=0;i<=9;i++))
do
echo "#$ -N bin$i" >L$i.sge
echo "#$ -S /bin/sh" >>L$i.sge
echo "#$ -cwd" >>L$i.sge
echo "/bin/echo Running on host: `hostname`.">>L$i.sge
echo "/bin/echo Starting on: `date`" >>L$i.sge
echo "cd /home/buldyrev/segregation/$1-$i" >>L$i.sge
echo "/bin/echo In directory: `pwd`" >>L$i.sge
echo "/home/buldyrev/segregation/f1.out < /home/buldyrev/segregation/start-f > ./Sout" >>L$i.sge 
qsub L$i.sge
done

