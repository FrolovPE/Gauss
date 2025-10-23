#!/bin/bash

if [ $# != 2 ]; then
	echo "Usage: supertest.sh <a.out> <r>"
else
    exe=$1
    r=$2
    
    if [ -f $exe ];then
    echo
    else
        echo "$exe doesn't exist"
        exit
    fi

    for((i=2; i<1000; i++)); do

        for((j=1; j<i; j++)); do

            for((k=1; k<4;k++)); do
                    echo "--------------------------------------------------------------------------"
                    echo "$exe $i $j $r $k"
                    $exe $i $j $r $k
                    echo ""
            done

        done

    done

fi