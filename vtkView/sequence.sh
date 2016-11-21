#!/bin/bash
# from small sequences of pictures makes one sequence
# ffmpeg video make preprocessing
#
# 14.8.2015
# Ondrej Vysocky - xvysoc01@stud.fit.vutbr.cz
# EPCC, University of Edinburgh

_help="	HELP:\n
	-n	number of sequences\n
	-i	input sequences name\n
	-o	output sequence name\n
	-f	image format\n
	-d	number of necessary digits for image index\n
	-h	help\n\n
if '-i name' then input sequences should be named:\n\t
	name1.0001.png, name15.0055.png...\n\t
	start indexing with 1"

# parse params
cnt=0
while getopts :n:i:o:f:d:h opt
do
	case "$opt" in
	n) cnt=$((cnt+1));nsq=$OPTARG ;;
	i) cnt=$((cnt+1));inp=$OPTARG ;;
	o) cnt=$((cnt+1));out=$OPTARG ;;
	f) cnt=$((cnt+1));typ=$OPTARG ;;
	d) cnt=$((cnt+1));dig=$OPTARG ;;
	h) echo -e $_help >&2; exit 1 ;;
	\?) echo -e $_help >&2; exit 1 ;;
	:) echo -e $_help >&2; exit 1 ;;
	esac
done

if [ $cnt -ne 5 ]
then
	echo -e  $_help >&2
	exit 1
fi

c=0
export i=1
while [ $c -ne $nsq ]
do
	for file in $inp$c.*
	do 
		x=$(printf "%0"$dig"d" $i)
		cp $file $out.$x.$typ
		i=$((i+1))
	done
	c=$((c+1))
done









