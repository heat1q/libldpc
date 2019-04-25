#!/bin/bash

codefile="$1"
simfile="$2"
mapfile="$3"
layerfile="$4"

if ! [ -e $codefile ] || [ "$codefile" = "" ]
then
    echo "Cannot open Codefile $codefile!"
    exit
fi

if ! [ -e $simfile ] || [ "$simfile" = "" ]
then
    echo "Cannot open simulation file $simfile!"
    exit
fi

if ! [ -e $mapfile ] || [ "$mapfile" = "" ]
then
    echo "Cannot open mapping file $mapfile!"
    exit
fi

if [ "$4" != "" ]
then
    if ! [ -e $layerfile ] 
    then
        echo "Cannot open layer file $layerfile!"
        exit
    fi
fi

#read sim name
while read -r -p "Enter name for simulation: " sim_name
do
    if [ -e $sim_name ] || [ "$sim_name" = "" ]
    then
        echo "Cannot use $sim_name. Simulation with this name already exists!"
    else
        break;
    fi
done

temp_dir="sim_$sim_name"

if ! [ -e $temp_dir ]
then
	mkdir $temp_dir
else
    if ! [ -z "$(ls -A $temp_dir)" ]
    then
        echo "Directory $temp_dir/ already exists!"
        exit
    fi
fi

#compile
sh make.sh $sim_name

cp "$simfile" "$temp_dir/sim.txt"
mv "$sim_name" "$temp_dir"

cd $temp_dir

#set the right name for the results
sed -i "1s/.*/name: res_$sim_name.txt/" sim.txt

#make executable
chmod +x $sim_name

#run
if [ "$4" != "" ]
then
    ./$sim_name -code "../$codefile" -sim sim.txt -map "../$mapfile" -layer "../$layerfile"
else
    ./$sim_name -code "../$codefile" -sim sim.txt -map "../$mapfile"
fi

#cleanup
rm "$sim_name"
rm "sim.txt"