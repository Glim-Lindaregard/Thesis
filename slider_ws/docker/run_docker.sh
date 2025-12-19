#!/bin/bash

#export DISPLAY=:0

xhost +local:root; 

IMG="ros2_slider_latest:nov2025"
NAME="ros2_slider_latest"

SKIP_FLAG="true"
PERSISTENT_FLAG="false"

if [[ "$1" == "--run-example" ]]; then
  SKIP_FLAG="false"
  shift
fi

if [[ "$1" == "--persistent" ]]; then
  PERSISTENT_FLAG="true"
  shift
fi

REMOVE_OPT="--rm"
if [[ "$PERSISTENT_FLAG" == "true" ]]; then
  REMOVE_OPT=""
fi

sudo docker run -it --rm $REMOVE_OPT \
  --name "$NAME" \
  --network=host \
  --privileged \
  --cap-add=NET_ADMIN \
  --device /dev/dri \
  -e SKIP_ENTRYPOINT=$SKIP_FLAG \
  -e DISPLAY=$DISPLAY \
  -e QT_QPA_PLATFORM=xcb \
  -e QT_X11_NO_MITSHM=1 \
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
  -v $HOME/.Xauthority:/root/.Xauthority:ro \
  -e XAUTHORITY=/root/.Xauthority \
  -v "$(eval echo ~$SUDO_USER)/slider_ws:/root/slider_ws" \
  -w /root/slider_ws \
  -v /lib/modules:/lib/modules:ro \
  "$IMG"
