# Slider ROS2 Workspace


ROS2 workspace for Slider zero-friction platform.


# Set up docker

Navigate to the docker sub-directory and build docker image with:

```bash 
docker build -t ros2_slider .
```

To run a docker container, run the *run_docker.sh* script:

```bash
./run_docker.sh --run-example
```

set --run-example option to run sample RL test scenario from entrypoint.sh 

In *run_docker.sh*, the workspace is mounted on to the container according to:

```bash
-v "$(eval echo ~$SUDO_USER)/slider_ws:/root/slider_ws" \
```


# Test RL Controller:


```bash
ros2 launch rl_controller go_2_point.launch.py \
  target_x:=1.2 \
  target_y:=1.1 \
  policy_file_path:=/root/slider_ws/src/slider/rl_controller/rl_controller/policies/policy.pth \
  record_bag:=false \
  launch_plotjuggler:=true \

```


```bash
ros2 launch rl_controller lissajous_trajectory_track.launch.py \
  lissa_omega:=0.07 lissa_A:=2 lissa_B:=2 lissa_a:=1 lissa_b:=1\
  record_bag:=false \
  launch_plotjuggler:=true \
  policy_file_path:=/root/slider_ws/src/slider/rl_controller/rl_controller/policies/policy.pth
```
