#!/bin/bash


bed=$1
isoforms=$2
output=$3

bed2gtf --bed $bed --isoforms $isoforms --output $output
