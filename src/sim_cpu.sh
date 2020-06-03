#!/bin/bash

codefile="$1"
simfile="$2"
mapfile="$3"
layerfile="$4"
defines="$5"

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

#read sim name
while read -r -p "Enter name for simulation: " sim_name
do
    if [ -e "../results/sim_"$sim_name"_cpu" ] || [ "$sim_name" = "" ]
    then
        echo "Cannot use $sim_name. Simulation with this name already exists!"
    else
        break;
    fi
done

sim_name=""$sim_name"_cpu"
temp_dir="sim_$sim_name"
if ! [ -e "../results" ]
then
    mkdir ../results
fi
mkdir "../results/$temp_dir"


#compile
echo "Compiling: sh make.sh $sim_name $defines"
sh make.sh $sim_name "$defines"

cp "$simfile" "../results/$temp_dir/sim.txt"
mv "$sim_name" "../results/$temp_dir"

cd "../results/$temp_dir"

#set the right name for the results
sed -i "1s/.*/name: res_$sim_name.txt/" sim.txt

#make executable
chmod +x $sim_name

#run#run
if [ "$4" != "" ]
then
    param="-code ../$codefile -sim sim.txt -map ../$mapfile -layer ../$layerfile"
else
    param="-code ../$codefile -sim sim.txt -map ../$mapfile"
fi

echo "Running: ./$sim_name $param"
./$sim_name $param

#cleanup
rm "$sim_name"
rm "sim.txt"
