You need CMAKE and Boost for bluild the program

> sudo apt-get install libboost-all-dev
> sudo apt-get -y install cmake

Create a directory build:

> mkdir build;

Go into that directory

> cd build;

You can choose to build a release version with:

> cmake ../ -DCMAKE_BUILD_TYPE=Release

or a debug version with

> cmake ../ -DCMAKE_BUILD_TYPE=Debug

And build everything with

> make

For run, go into that directory

> cd benchmarks

and run

> ./benchmark


To execute the benchmarks you will need python2.7 and python-tk

> sudo apt install python2.7 python-tk

Additionally, you need matplotlib in python2.7

> pip install --user matplotlib
