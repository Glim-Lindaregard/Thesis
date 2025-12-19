#!/bin/bash

isaaclabslider(){
    CACHE=/mimer/NOBACKUP/groups/ltu-rai-rl2024/nektaf/.isaac_cache
    DATA=/mimer/NOBACKUP/groups/ltu-rai-rl2024/nektaf/.isaac_data
    SLIDER_LIB_PATH=/cephyr/users/nektaf/Alvis/Desktop/mimer_ltu-rai-rl2024/nektaf/slider_isaac_lab/src
    mkdir -p "$CACHE" "$DATA"
    PUBLIC_IP=$(curl -s ifconfig.me)
    DISPLAY=$(env | grep DISPLAY=:)
    VGL_DISPLAY=$(env | grep VGL_DISPLAY=)

    apptainer shell --nv \
        --env $DISPLAY \
        --env $VGL_DISPLAY \
        --env "ACCEPT_EULA=Y" \
        --env LIVESTREAM=2 \
        --env YAML_C_EXT_DISABLED=1 \
        --env PYTHONNOUSERSITE=1 \
        --env OMNI_USER_DIR=/isaac-sim/kit/data \
        --env KIT_CACHE_DIR=/isaac-sim/kit/cache \
        --env KIT_USER_DIR=/isaac-sim/kit/data \
        --env "PYTHONPATH=/cephyr/users/nektaf/Alvis/.local/lib/python3.11/site-packages:${SLIDER_LIB_PATH}:${PYTHONPATH}" \
        --env "PATH=/cephyr/users/nektaf/Alvis/.local/bin:$PATH" \
        --bind "$CACHE":/isaac-sim/kit/cache \
        --bind "$DATA":/isaac-sim/kit/data \
        --bind /cephyr/users/nektaf/Alvis/Desktop/mimer_ltu-rai-rl2024/nektaf/isaac_overrides/manager_based_rl_env.py:/workspace/isaaclab/source/isaaclab/isaaclab/envs/manager_based_rl_env.py \
        /cephyr/users/nektaf/Alvis/Desktop/mimer_ltu-rai-rl2024/nektaf/containers/nektaf_isaaclab.sif
}

# --env $DISPLAY \
# --env $VGL_DISPLAY \
# DISPLAY=$(env | grep DISPLAY=:)
# VGL_DISPLAY=$(env | grep VGL_DISPLAY=)
