
rm smallpt_1 -rf

#g++ -O3 -fopenmp backup.cpp -o smallpt
#g++ -O3 -fopenmp smallpt.cpp -o smallpt -std=c++11 -ggdb

g++ -O0 -fopenmp smallpt.cpp -o smallpt_1 -std=c++11 -ggdb
#g++ -O0  smallpt.cpp -o smallpt -std=c++11 -ggdb

#g++ -O3 smallpt.cpp -o smallpt -std=c++11 -ggdb
#g++ -O0 smallpt.cpp -o smallpt -std=c++11 -ggdb
#g++ -O3 -fopenmp 00_smallpt.cpp -o smallpt 
if [ $? -ne 0 ]; then
	echo "error: compile src failed"
	exit -1;
fi

#exit 0;

rm image.ppm -rf
time ./smallpt_1 50
#time ./smallpt 100
#time ./smallpt 500
#time ./smallpt 5000
#time ./smallpt 50000

diff -q image_5K.ppm image.ppm



# 1. 插入三角形
# 2. 插入其他曲面
# 3. 用自己的函数渲染

# 3. 插入蝴蝶？
# 4. 插入桃花花瓣？
