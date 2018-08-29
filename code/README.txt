<<<<<<< HEAD
Unzip the project, open a terminal and change your working directory to the unzipped folder
=======
You need CMAKE and Boost for bluild the program

> sudo apt-get install libboost-all-dev
> sudo apt-get -y install cmake
>>>>>>> e9674d98bd2ff8f46a2fcacc0ef26c5425b0269a

Create a directory build:

> mkdir build;

Go into that directory

> cd build;

You can choose to build a release version with:

> cmake ../code -DCMAKE_BUILD_TYPE=Release

or a debug version with

> cmake ../code -DCMAKE_BUILD_TYPE=Debug

And build everything with

> make

<<<<<<< HEAD
To execute the tests go to the /build/benchmarks directory

> cd build/benchmarks

To execute the tests you can use

> ./benchmark -t RootFinders -r detailed

To execute the test that shows tha plots you can use

> ./benchmark -t RootFindersPlotted

RootFindersPlotted will show plots showing the amount of test function calls
for each epsilon, for ech of the root finding methods. It shows one plot after
the other for each of the test functions 


**********************************************************************************8
******* Dependencies **********************************************************
=======
For run, go into that directory

> cd benchmarks

and run

> ./benchmark

>>>>>>> e9674d98bd2ff8f46a2fcacc0ef26c5425b0269a

To execute the benchmarks you will need python2.7 and python-tk

> sudo apt install python2.7 python-tk

Additionally, you need matplotlib in python2.7

> pip install --user matplotlib


