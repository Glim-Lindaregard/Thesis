#!/usr/bin/env bash
set -e

WS=~/slider_ws
BAG_DIR="$WS/rosbags"
BAG_NAME="ftc_test"

mkdir -p "$BAG_DIR"
rm -rf "$BAG_DIR/$BAG_NAME"

cd "$WS"
colcon build
source "$WS/install/setup.bash"

cleanup() {
  echo "Stopping thrust..."
  ros2 topic pub --once /eight_thrust_pulse std_msgs/msg/UInt8MultiArray \
  "{data: [0,0,0,0,0,0,0,0]}"

  echo "Stopping ROS processes..."
  kill ${PID_BAG:-} ${PID_FTC:-} ${PID_TARGET:-} ${PID_GAZEBO:-} 2>/dev/null || true

  # give bag a moment to flush
  wait ${PID_BAG:-} 2>/dev/null || true

  echo "Plotting latest bag..."
  python3 "$WS/src/slider/slider_ftc_control/slider_ftc_control/viz/plot_states.py" \
    "$BAG_DIR/$BAG_NAME/${BAG_NAME}_0.db3" || true
}
trap cleanup INT TERM EXIT

echo "Stopping thrust..."
ros2 topic pub --once /eight_thrust_pulse std_msgs/msg/UInt8MultiArray \
  "{data: [0,0,0,0,0,0,0,0]}"

#echo "Starting Gazebo..."
#ros2 launch slider_gazebo slider_launch.py gui:=true &
#PID_GAZEBO=$!

sleep 5

echo "Starting Vicon..."
ros2 launch mocap_pose_to_odom get_odom_vicon.launch.py &

sleep 5

echo "Starting target trajectory..."
ros2 launch target_trajectories target_point.launch.py \
  target_x:=0.0 target_y:=0.0 target_theta:=0.0&
PID_TARGET=$!
sleep 2

echo "Starting FTC node..."
ros2 run slider_ftc_control ftc_node --ros-args -p pwm_resolution:=64 &
PID_FTC=$!

sleep 0.2

echo "Starting rosbag recording to $BAG_DIR/$BAG_NAME ..."
ros2 bag record -o "$BAG_DIR/$BAG_NAME" \
  /odom \
  /target_point \
  /eight_thrust_pulse \
  /thruster_state &
PID_BAG=$!

echo "All processes started."
echo "PIDs: Gazebo=$PID_GAZEBO Target=$PID_TARGET FTC=$PID_FTC Bag=$PID_BAG"
echo "Press Ctrl+C to stop everything."

wait
