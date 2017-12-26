echo 1
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K1 -lambdaweight 1 -ic 300 -sic 1000 >logK1.txt 2>&1
echo 2
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K2 -lambdaweight 2 -ic 300 -sic 1000 >logK2.txt 2>&1
echo 4
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K4 -lambdaweight 4 -ic 300 -sic 1000 >logK4.txt 2>&1
echo 8
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K8 -lambdaweight 8 -ic 300 -sic 1000 >logK8.txt 2>&1
echo 16
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K16 -lambdaweight 16 -ic 300 -sic 1000 >logK16.txt 2>&1
echo 24
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K24 -lambdaweight 24 -ic 300 -sic 1000 >logK24.txt 2>&1
echo 32
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K32 -lambdaweight 32 -ic 300 -sic 1000 >logK32.txt 2>&1
echo 48
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K48 -lambdaweight 48 -ic 300 -sic 1000 >logK48.txt 2>&1
echo 64
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K64 -lambdaweight 64 -ic 300 -sic 1000 >logK64.txt 2>&1
echo 128
./main -i ../../Data/weibo/ -o ../../Result/weibo/CTCD/K128 -lambdaweight 128 -ic 300 -sic 1000 >logK128.txt 2>&1


