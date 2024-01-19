# AgentSeal 2.0 - Grey Seals and Harbor Seals

## To get started

1. The Model is coded in NetLogo 6.4.0 available to download for free [here](https://ccl.northweatern.edu/netlogo/). You can download for Windows, Mac, or Linus OS. See the [NetLogo User Manual](https://ccl.northwestern.edu/netlogo/docs/) for more information about the software.


2. Clone the [AgentSeal_with_GreySeals]() GitHub repository onto your local system.
    
    Enter the below code into your computer terminal.

    `git clone https://github.com/KaraWatts/AgentSeal_with_GreySeals.git`

3. AgentSeal2.0_with_GreySeals.nlogo is the netlogo file containing GUI and the code so once you download NetLogo, double-clock this file and it should load in NetLogo.

## Running the model

1. Use the sliders to set the n_of_gseals (# of grey seals) and n_of_hseals (# of harbor seals) to the amounts you want to model. For more precision you can left click on each input block, select edit, and imput the exact value under 'Value'.

2. Set the monthsOfSimulation to the number of months you wish to simulate.

3. Click 'setups' to render the map. This will include all haulout sites and grey and harbor seals distributed across haulouts.

4. Click ![go button](/AgentSeal_with_GreySeals/go.png) to progress the model by individual steps OR click ![go repeating](/AgentSeal_with_GreySeals/go_repeating.png) to run continuously until the end of the model duration. (You can end a model run early by clicking ![go repeating](/AgentSeal_with_GreySeals/go_repeating.png) again)

## Folder Breakdown

### OutPut Folder

This will be where all model output files are saved. Make sure to move any output files into new folders before running a new model session or it will save over existing files.

### Input Folder

Input folder contains files necessary in two Procedures: `set_landscape` and `load_haulouts`.
- `set_landscape`:
    - *GrecianEtAl2018.asc* – raster with HSI values of each at-sea patch
    - *land_distance_EastCoast_noBay_large.asc* – raster file showing distance from each at-sea patch to the coast
    - *land_EastCoast_5km_NoBay_large.asc* and *land_EastCoast_25km_NoBay_large.asc and* – 5x5 km and 25x25 km grid ids. Necessary for memory procedures
    - *ho_distance_id_GS_X.asc* – 9 files each showing distance of each patch to a given haul-out site X
    - *Maxhsi_EastCoast_25km_NoBay_large.asc* – maximum initial value of HSI for each 25x25km grid
-	`load_haulouts`:
    - *HaulOutPoints_EastCoast_GreySeal.txt* – location of 9 haul-outs sites with state variables of these sites

## Interface

![AgentSeal Interface](/AgentSeal_with_GreySeals/Interface_Breakdown.jpeg)



