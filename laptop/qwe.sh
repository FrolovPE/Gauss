#!/bin/bash

read -p "asd: " f

echo "$f"


if [ $f != "cc" ]; then
	echo "Usage <arg>"
else
	echo nice
fi
