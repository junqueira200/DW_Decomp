#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo "Error: No argument provided!"
    echo "Usage: $0 <folder>"
    exit 1
fi

cp -r build2 $1
mkdir $1/CSV
mkdir $1/Plot
mkdir $1/Plot/Heuristic
mkdir $1/Plot/Model
EXEC=$1/build2/2L_SDVRP_MILP_Main/2L_SDVRP_MILP_Main
#echo $EXECC
CURRENT_DATE=$(date)
echo "# Created at: $CURRENT_DATE" > "$1/oroloc3D.csv"
echo "Inst; Seed; TruckId; NumberOfItems; ModelStatus; ModelTime; ExtremePoitsHeurStatus; ExtremePoitsHeurTime; GA_Status; GA_Time" >> "$1/oroloc3D.csv"

#CONFIG="--compactness 0.8 --axleWeights 0.37 --support 0.12 --lifo 0.53 --itemsNotPacked 0.62 --loadBalancing 0.92 --gaDiscard 0.15 --gaElite 0.18 --gaIslands 4 --gaMigrInterval 72 --gaMigtRate 0.02 --gaInherit 0.72"

#CONFIG="--compactness 0.02 --axleWeights 0.23 --support 0.57 --lifo 0.02 --itemsNotPacked 1.31 --loadBalancing 0.02 --gaDiscard 0.12 --gaElite 0.2 --gaIslands 1"

CONFIG="--compactness 0.06 --axleWeights 1.19 --support 7.6 --lifo 2.4 --itemsNotPacked 100.0 --loadBalancing 0.6 --gaDiscard 0.3 --gaElite 0.26 --gaIslands 4 --gaMigrInterval 38 --gaMigtRate 0.02 --gaInherit 0.69 --gaSoftGroup 1.52"

./$EXEC -i TransforData/Instances/S10_240T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S10_240T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S20_288T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S20_288T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S30_288T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S30_288T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S40_336T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S40_336T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S50_360T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S50_360T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S60_360T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S60_360T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S70_384T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S70_384T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S80_408T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S80_408T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S90_432T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S90_432T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S96_456T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S96_456T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
./$EXEC -i TransforData/Instances/S99_480T.txt --solOroloc3D_2 TransforData/TrasnforSolution/S99_480T_sol_items.csv --oroloc3D_2 --strOutputFolder $1 $CONFIG
