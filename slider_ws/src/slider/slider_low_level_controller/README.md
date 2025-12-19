# slider_low_level_controller


Converts continuous body frame thrust and torque commands to discrete PWM thrust pulses 

Contents: 
- thrust_to_pwm.py -> converts (3,) force vector to (8,) thrust vector for each on-board thruster
- write_to_arduino(for real robot deployment) -> subscribes to (8,) thrust vector and writes data to arduino board with specified baud rate 

Use launch files to launch the ros2 nodes 



