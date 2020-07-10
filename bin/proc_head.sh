#!/usr/bin/env bash

awk '{if (NR%4==1) {if ($0~/p3-/) a=gensub(/[^p]* (p3-[_0-9A-Z-]+) .*/,"\\1",$0); else a="NONE"; print "@"a,$1,$2} else print}' $1
