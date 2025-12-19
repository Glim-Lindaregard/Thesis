#!/bin/bash

set -e

echo "Starting..."

if [[ "${SKIP_ENTRYPOINT:-false}" == "true" ]]; then
  echo "Skipping entrypoint setup."
  exec bash -i
fi


source /opt/ros/humble/setup.bash

colcon build

source install/setup.bash

ros2 launch rl_controller go_2_point.launch.py \
  target_x:=5.2 \
  target_y:=1.1 \
  policy_file_path:=/root/slider_ws/src/slider/rl_controller/rl_controller/policies/policy.pth \
  frame_id:=odom \
  publish_rate:=0

echo "Done"

if [[ $# -gt 0 ]]; then
  # run user command in the already-prepared env
  exec "$@"
else
  # keep the same env; just open an interactive shell
  exec bash -i
fi

