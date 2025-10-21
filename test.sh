#!/bin/bash


if [ $# != 1 ]; then
	echo "Usage: test.sh <dir>"
else
	 for i in $1*; do
	     ii="${i##*/}"
	     filename="${i%.*}"
	     filename="${filename##*/}"
	     ext="${i##*.}"
#	     echo "$ext"
#	     echo $i
		if [ $ext != "sh" ] && [ $ii != "3" ] && [ $ext != "." ] && [ $ext != "zip" ]; then
		         if [ $filename == "c" ] || [ $filename == "d" ] || [ $filename == "e" ]; then
				 echo $i
		                 echo "./a.out 6 2 6 0 $ii"
		                 ./a.out 6 2 6 0 $i
		                 echo ""
		         else
		             echo ""
		             echo $i
			     echo "./a.out 4 2 4 0 $ii"
		             ./a.out 4 2 4 0 $i
		         fi
		fi
	 done
fi
