#!/bin/bash
cd build
cmake D:/Study/Codes/AnglesC -G "MinGW Makefiles"
mingw32-make main
./main.exe