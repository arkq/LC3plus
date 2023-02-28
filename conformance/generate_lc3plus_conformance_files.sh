#!/bin/bash

case $# in
    1) tag=$1;;
    *) echo usage: $0 tag; exit 1;;
esac

timestamp=$(date +"%Y%m%d_%H%M%S")

config_names="LC3plus_HighResolution_Streaming LC3plus_HighResolution_fallback_rates LC3plus_DECT_AAP_HighResolution LC3plus_DECT_AAP_service LC3plus_DECT_voice_service"
generated_files=""
config_files=""
for config in $config_names; do
    generated_files="$generated_files test_files_${config}/ ${config}*.log"
    config_files="$config_files ${config}.cfg"
done
echo removing files from previous run
rm -rf $generated_files

ret_val=0
for config in $config_names; do
    echo python3 lc3plus_conformance.py -wav_only ${config}.cfg
    python3 lc3plus_conformance.py -wav_only ${config}.cfg
    ret_val=$(($ret_val + $?))
done
if [ $ret_val == 0 ]; then
    echo  
    echo file generation with lc3plus_conformance.py successfull!
else
    echo  
    echo file generation with lc3plus_conformance.py failed! see logs for errors!
    exit 1
fi

zip lc3plus_conformance_files${timestamp}_${tag}.zip -r $generated_files $config_files
ret_val=$(($ret_val + $?))

if [ $ret_val == 0 ]; then
    echo  
    echo Process complete! 
    echo wave files zipped to lc3plus_conformance_files${timestamp}_${tag}.zip
else
    echo  
    echo zipping failed!
    exit 1
fi