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

        DeclareLaunchArgument('lissa_A', default_value='2.0', description='Amplitude A (float)'),
        DeclareLaunchArgument('lissa_B', default_value='2.0', description='Amplitude B (float)'),
        DeclareLaunchArgument('lissa_a', default_value='1', description='Frequency a (int)'),
        DeclareLaunchArgument('lissa_b', default_value='2', description='Frequency b (int)'),
        DeclareLaunchArgument('lissa_delta', default_value='1.0', description='Phase delta (float)'),
        DeclareLaunchArgument('lissa_omega', default_value='0.05', description='Angular frequency (float)'),

        Node(
            package='target_trajectories',
            executable='lissajous_trajectory',
            name='lissajous_trajectory',
            parameters=[{
                'lissa_A': LaunchConfiguration('lissa_A'),
                'lissa_B': LaunchConfiguration('lissa_B'),
                'lissa_a': LaunchConfiguration('lissa_a'),
                'lissa_b': LaunchConfiguration('lissa_b'),
                'lissa_delta': LaunchConfiguration('lissa_delta'),
                'lissa_omega': LaunchConfiguration('lissa_omega'),
                'frame_id': 'odom',
                'publish_rate': 10.0,
            }],
            output='screen',
        ),
    ])
