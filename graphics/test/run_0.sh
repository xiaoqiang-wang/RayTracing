
rm smallpt_0 -rf

#g++ -O3 -fopenmp backup.cpp -o smallpt
#g++ -O3 -fopenmp smallpt.cpp -o smallpt -std=c++11 -ggdb

#g++ -O0 -fopenmp smallpt.cpp -o smallpt -std=c++11 -ggdb
g++ -O0  smallpt.cpp -o smallpt_0 -std=c++11 -ggdb

#g++ -O3 smallpt.cpp -o smallpt -std=c++11 -ggdb
#g++ -O0 smallpt.cpp -o smallpt -std=c++11 -ggdb
#g++ -O3 -fopenmp 00_smallpt.cpp -o smallpt 
if [ $? -ne 0 ]; then
	echo "error: compile src failed"
	exit -1;
else
	echo "OK: compile src success."
	#exit 0;
fi

#exit 0;

export DEBUG_RT=1

unset ENALBE_VALUE_CHECK
export ENALBE_VALUE_CHECK=1


rm rt_dump* -rf

rm image.ppm -rf
time ./smallpt_0 50
#time ./smallpt 100
#time ./smallpt 500
#time ./smallpt 5000
#time ./smallpt 50000

grep -r "radiance_0"  rt_dump.txt >rt_dump_0.txt
grep -r "radiance_1"  rt_dump.txt >rt_dump_1.txt
vimdiff rt_dump_0.txt rt_dump_1.txt

#diff -q image_5K.ppm image.ppm



# 1. 插入三角形
# 2. 插入其他曲面
# 3. 用自己的函数渲染

# 3. 插入蝴蝶？
# 4. 插入桃花花瓣？
