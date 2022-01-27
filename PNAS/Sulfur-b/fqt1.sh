for((j=0;j<=9;j+=1))
do
for((i=0;i<=9;i++))
do
echo "#$ -N fq$j-$i" >L$j-$i.sge
echo "#$ -S /bin/sh" >>L$j-$i.sge
echo "#$ -cwd" >>L$j-$i.sge
echo "/bin/echo Running on host: `hostname`.">>L$j-$i.sge
echo "/bin/echo Starting on: `date`" >>L$j-$i.sge
echo "cd /home/buldyrev/segregation/$1$2.$j-$i" >>L$j-$i.sge
echo "/bin/echo In directory: `pwd`" >>L$j-$i.sge
echo "/home/buldyrev/segregation/f1.out < /home/buldyrev/segregation/start-f > ./Sout" >>L$j-$i.sge 
qsub L$j-$i.sge
done
done

