#!/usr/bin/env bash
set -e

########### Script to build relevant packages, source files, start nodes, start recording bags, and plot bags ############

# ------ YOU CAN ALSO SET SET POINT HERE ---

# --- config ---
WS=~/slider_ws
BAG_DIR="$WS/rosbags"
BAG_NAME="ftc_test"

mkdir -p "$BAG_DIR"
rm -rf "$BAG_DIR/$BAG_NAME"

cd "$WS"
colcon build --packages-select slider_gazebo slider_ftc_control target_trajectories

# Source workspace
source "$WS/install/setup.bash"


echo "Starting Gazebo + slider..."
ros2 launch slider_gazebo slider_launch.py &
PID_GAZEBO=$!

# Give Gazebo some time to come up
sleep 5

echo "Starting target trajectory..."

################## SET SETPOINT HERE ####################
ros2 launch target_trajectories target_point.launch.py \
  target_x:=1.0 target_y:=0.4  target_theta:=0.4 &
PID_TARGET=$!
#########################################################
sleep 2

echo "Starting FTC node..."
ros2 run slider_ftc_control ftc_node --ros-args -p failure_time_s:=2.0 &
PID_FTC=$!

sleep 0.1

echo "Starting rosbag recording to $BAG_DIR/$BAG_NAME ..."
ros2 bag record -o "$BAG_DIR/$BAG_NAME" \
    /odom \
    /target_point \
    /eight_thrust_pulse \
    /thruster_state &
PID_BAG=$!

echo "All processes started."
echo "PIDs: Gazebo=$PID_GAZEBO  Target=$PID_TARGET  FTC=$PID_FTC  Markers=$PID_MARKERS  Bag=$PID_BAG"
echo "Press Ctrl+C to stop everything."

cleanup() {
    echo "Stopping ROS processes..."
    # Try to kill gracefully; ignore errors if some already died
    kill $PID_BAG $PID_MARKERS $PID_FTC $PID_TARGET $PID_GAZEBO 2>/dev/null || true

    # ensure bag file is fully flushed
    wait $PID_BAG 2>/dev/null || true

    echo "Plotting latest bag..."
    python3 "$WS/src/slider/slider_ftc_control/slider_ftc_control/viz/plot_states.py" \
        "$BAG_DIR/$BAG_NAME/${BAG_NAME}_0.db3" || true
}

trap cleanup INT TERM

# wait for background processes until Ctrl+C
wait
