from launch import LaunchDescription
from launch_ros.actions import Node
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch_ros.substitutions import FindPackageShare
from launch.actions import IncludeLaunchDescription, DeclareLaunchArgument, ExecuteProcess
from launch.substitutions import LaunchConfiguration
from launch.conditions import IfCondition
import os
from datetime import datetime
from ament_index_python.packages import get_package_share_directory


high_level_path = os.path.join(os.getcwd(), 'rosbags')
os.makedirs(high_level_path, exist_ok=True)

bag_name = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
bag_path = os.path.join(high_level_path, bag_name)

def generate_launch_description():
    return LaunchDescription([

        DeclareLaunchArgument('target_x', default_value='0.0', description=''),
        DeclareLaunchArgument('target_y', default_value='0.0', description=''),
        DeclareLaunchArgument('target_theta', default_value='0.0', description=''),

        Node(
            package='target_trajectories',  
            executable='target_point',
            name='target_point',
            parameters=[{
                'frame_id': 'odom',
                'publish_rate': 10.0,
                'target_x': LaunchConfiguration('target_x'),
                'target_y': LaunchConfiguration('target_y'),
                'target_theta': LaunchConfiguration('target_theta'),
            }],
        )
        
    ])
