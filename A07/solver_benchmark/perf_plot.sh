for i in $(ls)
do 
numcores=$(echo $i | awk -F_ '{print $1}')
cd "$i"
perf=$(cat resfile | grep "Performance total:" resfile | tail -n 1 | awk '{print $(NF-1)}')
echo "$numcores $perf"
cd ..
done
