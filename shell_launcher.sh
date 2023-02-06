#!/bin/bash
#
#
#SIMULATION PARAMETERS:
#
#Name of the simulation
nameFile="_EXAMPLE"
#Number of cores used (for matrix inversion, max 32)
nCores=32
#Number of threads used (to fill in the matrices)
nThreads=64
#
#
#RAM/CPU MONITOR:
#
#Options to monitor in real time the RAM and CPU used by the simulation
#Set the following variable to the values: 1 to monitor the resources, or 0 to do not
yesORnoMonitor=1
maxHoursRamMonitor=240
#
#
#
#LAUNCHING THE SIMULATION:
#
echo "Launching Julia files"
#
#nohup /usr/bin/time -v 
nohup julia -p $nCores -t $nThreads "Launcher.jl" $nameFile > "Data_Output/out$nameFile"_p"$nCores"_t"$nThreads.out" &
codeSimulation=$!
#
if [ $yesORnoMonitor == 1 ]
then
    nohup julia -p 1 "Resources_Monitor/Resources_Monitor.jl" $nameFile $codeSimulation $maxHoursRamMonitor > "Resources_Monitor/res_out$nameFile"_p"$nCores"_t"$nThreads.out" &
    echo "Monitoring the resources"
fi
#
echo "Launching completed"