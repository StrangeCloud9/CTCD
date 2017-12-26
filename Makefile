exe = main
exe_poi = main_poi
cc = g++
objects = ./dataset.o ./main.o ./model.o ./utils.o
objects_poi = ./dataset_poi.o ./main_poi.o ./model_poi.o ./utils_poi.o
std = c++11

all : $(exe) $(exe_poi)

$(exe) : $(objects)
			$(cc) -O2 -fopenmp -std=$(std) -o $(exe) $(objects)
$(exe_poi) : $(objects_poi)
			$(cc) -O2 -fopenmp -std=$(std) -o $(exe_poi) $(objects_poi)

./dataset.o : ./dataset.cc dataset.h
	$(cc) -std=$(std) -O2 -fopenmp -c ./dataset.cc -o ./dataset.o
./main.o : ./main.cc utils.h dataset.h main.h model.h
	$(cc) -std=$(std) -O2 -fopenmp -c ./main.cc -o ./main.o
./model.o : ./model.cc utils.h model.h dataset.h
	$(cc) -std=$(std) -O2 -fopenmp -c ./model.cc -o ./model.o
./utils.o : ./utils.cc utils.h
	$(cc) -std=$(std) -O2 -fopenmp -c ./utils.cc -o ./utils.o

./dataset_poi.o : ./dataset.cc dataset.h
	$(cc) -std=$(std) -O2 -DPOISSON -fopenmp -c ./dataset.cc -o ./dataset_poi.o
./main_poi.o : ./main.cc utils.h dataset.h main.h model.h
	$(cc) -std=$(std) -O2 -DPOISSON -fopenmp -c ./main.cc -o ./main_poi.o
./model_poi.o : ./model.cc utils.h model.h dataset.h
	$(cc) -std=$(std) -O2 -DPOISSON -fopenmp -c ./model.cc -o ./model_poi.o
./utils_poi.o : ./utils.cc utils.h
	$(cc) -std=$(std) -O2 -DPOISSON -fopenmp -c ./utils.cc -o ./utils_poi.o

clean :
		rm $(objects) $(exe) $(objects_poi) $(exe_poi)
