import sys
if sys.prefix == '/usr':
    sys.real_prefix = sys.prefix
    sys.prefix = sys.exec_prefix = '/home/glim/Thesis/Code/Gazebo/install/slider_thruster_controller'
