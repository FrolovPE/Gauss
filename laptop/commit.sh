#!/bin/bash

export asd="$(date -R | head -c 26)"

git add main.cpp lib.cpp lib.h makefile
git commit -m $asd
echo $asd > lastcom
	 
