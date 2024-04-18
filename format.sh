#!/bin/bash

# SLURM options:

#SBATCH --job-name=format                         # Job name
#SBATCH --output=format.log                # Standard output and error log
#SBATCH --partition=htc                              # Partition choice
#SBATCH --ntasks=1                                   # Maximum number of parallel processes
#SBATCH --mem=200G                                    # Amount of memory required
#SBATCH --time=24:00:00                              # 7 days by default on htc partition


###################################
# Set up environment

source ~/.bashrc
source $HOME/Envs/setupenv_WCSim-1p12p1.sh
# source $HOME/Envs/setupenv_WCSim-sol..sh


###################################
# Var declaration

dir=/sps/t2k/lperisse/Soft/Preprocess_WCSim
darkrate=4.2                # PMT darkrate in kHz
minimum_hit_sig=1           # Minimum number of hits of signal in a trigger window to consider it a signal windowed event
max_event_bg=0              # Maximum number of background event to extract
nhit_per_200ns=10            # SK-SLE: 25
eventtype=1
verbose=1

#HK
tank_rad=3240.0               #HK:3240.0,    SK:1684.0750, WCTE:172.05
tank_halfz=3287.55            #HK:3287.55,   SK:1810.0,    WCTE:136.95
nb_PMT_ID=19547               # Number of PMTs in the ID of the detector,  SK:11129,  HK:19547

#SK
# tank_rad=1684.0750          #HK:3240.0,    SK:1684.0750, WCTE:172.05
# tank_halfz=1810.00          #HK:3287.55,   SK:1810.0,    WCTE:136.95
# nb_PMT_ID=11129             # Number of PMTs in the ID of the detector,  SK:11129,  HK:19547


sig_file=/sps/t2k/lperisse/Soft/wcsim/results/wcsim112_UnifVtx_electron_HK.root
# sig_file=/sps/t2k/lperisse/Soft/wcsim/results/wcsim_solar_HK_darkrateOFF_anytime.root
# sig_file=/sps/t2k/lperisse/Soft/wcsim/results/wcsim_solar_HK_darkrateOFF_anytime.root
# sig_file=/sps/t2k/lperisse/Soft/wcsim/results/wcsim_solar_HK_darkrateOFF_daytime.root
# sig_file=/sps/t2k/lperisse/Soft/wcsim/results/wcsim_solar_HK_darkrateOFF_nighttime.root


noise_file=/sps/t2k/lperisse/Soft/wcsim/results/wcsim112_DarkNoise_HK.root
# noise_file=/sps/t2k/lperisse/Soft/wcsim/results/wcte112_darknoise.root


###################################
# Commands to be submitted

echo "---------------------------------------------------------------------------------------"
echo ""

#Run AnalyzeWCSim
cd $dir
./format_data_from_WCSim  -n $noise_file  -s $sig_file  -e $eventtype  -r 2000  -t $tank_rad  -z $tank_halfz  -p $nb_PMT_ID  -d $darkrate  -w $minimum_hit_sig  -h $nhit_per_200ns  -v $verbose  -m $max_event_bg
# ./format_data_from_WCSim                -s $sig_file  -e $eventtype  -r 0     -t $tank_rad  -z $tank_halfz  -p $nb_PMT_ID  -d $darkrate  -w $minimum_hit_sig  -h $nhit_per_200ns  -v $verbose

echo ""
echo "---------------------------------------------------------------------------------------"
echo "Job done"


