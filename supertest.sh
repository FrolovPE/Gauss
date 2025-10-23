#!/bin/bash

if [ $# != 3 ]; then
	echo "Usage: supertest.sh <a.out> <r> <n(до какого n делать тесты)>"
else
    exe=$1
    r=$2
    n=$3
    
    if [ -f $exe ];then
    echo
    else
        echo "$exe doesn't exist"
        exit
    fi

    for((i=2; i<=$n; i++)); do

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