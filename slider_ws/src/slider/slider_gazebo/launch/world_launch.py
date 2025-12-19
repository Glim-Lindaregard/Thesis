import os

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.actions import IncludeLaunchDescription
from launch.conditions import IfCondition
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import LaunchConfiguration
from launch.actions import ExecuteProcess


def generate_launch_description():
    world_file = "empty_world.world"

    slider_package = "slider_gazebo"

    world_path = os.path.join(get_package_share_directory(slider_package), "worlds", world_file)

    return LaunchDescription([        
        ExecuteProcess(
            cmd=['gazebo', '--verbose', '-s', 'libgazebo_ros_factory.so', world_path],
            output='screen'),
    ])