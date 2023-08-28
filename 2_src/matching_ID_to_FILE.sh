#!/usr/bin/env bash
######################################
#		08/08/2023
#		By Elyna Bouchereau
######################################
for line in $(cut -f1 -d'|' $1); do grep $line $2; done