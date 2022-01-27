for((j=0;j<=9;j+=1))
do
for((i=0;i<=9;i++))
do
cd $1-$2.$j-$i
nice ../f1.out <../start-f >junk-out
cd ..
done
done

