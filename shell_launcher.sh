#!/bin/bash
#
#
#SIMULATION PARAMETERS:
#
#Name of the simulation
nameFile="_default"
#Number of cores used (for matrix inversion)
nCores=32
#Number of threadsused (for filling in the matrices)
nThreads=$nCores
#
#
#RAM/CPU MONITOR:
#
#Options to monitor in real time the RAM and CPU used by the simulation
#Set the following to "y" to monitor, or "n" to do not
yesORnoMonitor = "y"
maxHoursRamMonitor = 240
#
#
#
#LAUNCHING THE SIMULATION:
#
echo "Launching Julia files"
#
#nohup /usr/bin/time -v 
nohup julia -p $nCores -t $nThreads "Launcher.jl" $nameFile > "Data_Output/out$nameFile"_"$nCores"_"$nThreads.out" &
codeSimulation=$!
#
if [yesORnoMonitor == "y"]
    nohup julia -p 1 "RAM_monitor/RAM_monitor.jl" $nameFile $codeSimulation $maxHoursRamMonitor > "ram_out.out" &
fi
#
echo "Launching completed"