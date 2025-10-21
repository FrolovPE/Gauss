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
				 		
						for((k = 1; k <= 6 ; k++)); do
						echo "--------------------------------------------------------------------------"
		                 echo "./a.out 6 $k 6 0 $ii"
		                 ./a.out 6 $k 6 0 $i
		                 echo ""
						 done
		         else
		             echo ""
		             echo $i
						for((k = 1; k <= 4 ; k++)); do
						echo "--------------------------------------------------------------------------"
			     		echo "./a.out 4 $k 4 0 $ii"
		                ./a.out 4 $k 4 0 $i
						done
		         fi
		fi
		echo
	 done
fi
