#!/bin/bash

export asd="$(date -R | head -c 26)"
echo $asd > lastcom
export aaa="$(cat lastcom)"

git add main.cpp lib.cpp lib.h makefile test.sh
git commit -m $aaa
