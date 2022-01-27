for((i=0;i<=9;i++))
do
for((j=0;j<=9;j+=1))
do
echo "bcp_data 3 7 1 $1-$2.$i-$j-dc" > start-dc
cd $1-$2.$i-$j
../dc.out <../start-dc
cd ..
done
done
