#!/bin/bash
cd build
cmake C:/Users/redmi/Desktop/AnglesC -G "MinGW Makefiles"
mingw32-make main
./main.exe