#!/bin/bash

AGE_1=`stat -c %Y "$1"`
AGE_2=`stat -c %Y "$2"`

#echo ${AGE_1}
#echo ${AGE_2}

if [ "${AGE_1}" -gt "${AGE_2}" ]; then
	echo "No Changes"
else
	echo "Changes"
fi
