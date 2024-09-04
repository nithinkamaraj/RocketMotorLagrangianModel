#!/bin/bash

find -maxdepth 1 '(' -name "0.*" -o -name "1.*" -o -name "*e-*" ')' -print > tobeDeleted
for (( i=0 ; i<$1 ; i++ )); 
do
	cd processor$i/
	xargs rm -rf < ../tobeDeleted
	cd ..
	echo processor$i
done
# xargs rm -rf < tobeDeleted
