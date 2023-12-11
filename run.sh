#!/bin/bash

gcc -o simple_filter.exe simple_filter.c -lm

M=1024  

./simple_filter $M hL.txt hR.txt YL.txt YR.txt blue_giant_fragment.wav output.wav

python3 show_data.py
