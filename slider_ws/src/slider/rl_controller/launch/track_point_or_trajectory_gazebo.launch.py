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
    pkg_share = get_package_share_directory('rl_controller')
    plotjuggler_config = os.path.join(
        pkg_share, 'configs', 'plot_juggler_config_rl_target_tracking.xml'
    )

    return LaunchDescription([
        IncludeLaunchDescription(
            PythonLaunchDescriptionSource(
                [FindPackageShare('slider_gazebo'), '/launch/slider_launch.py']
            ),
        ),

        DeclareLaunchArgument('policy_file_path', description='policy_file_path'),
        DeclareLaunchArgument('record_bag', default_value='false', description='Set true to record rosbag'),
        DeclareLaunchArgument('bag_name', default_value='test', description='Output bag name'),
        DeclareLaunchArgument('launch_plotjuggler', default_value='false', description='Launch PlotJuggler'),

        ExecuteProcess(
            condition=IfCondition(LaunchConfiguration('record_bag')),
            cmd=[
                'ros2', 'bag', 'record',
                '--storage', 'sqlite3',
                '-o', bag_path,
                '/odom',
                '/target_point',
                '/action_bounds',
                '/thrust_cmd',
            ],
            output='screen',
        ),

        IncludeLaunchDescription(
            PythonLaunchDescriptionSource(
                [FindPackageShare('slider_low_level_controller'), '/launch/thrust_to_pwm.launch.py']
            ),
        ),

        Node(
            package='rl_controller',
            executable='policy',
            name='policy',
            arguments=['--ros-args', '--log-level', 'policy:=DEBUG'],
            parameters=[{
                'policy_file_path': LaunchConfiguration('policy_file_path'),
            }],
            output='screen',
        ),

        Node(
            package='rl_controller',
            executable='action_bounds',
            name='action_bounds',
            output='screen',
        ),
        
        
        IncludeLaunchDescription(
            PythonLaunchDescriptionSource(
                [FindPackageShare('target_trajectories'), '/launch/lissajous_trajectory.launch.py']
            ),
        ),
        
        # IncludeLaunchDescription(
        #     PythonLaunchDescriptionSource(
        #         [FindPackageShare('target_trajectories'), '/launch/target_point.launch.py']
        #     ),
        # ),

        ExecuteProcess(
            condition=IfCondition(LaunchConfiguration('launch_plotjuggler')),
            cmd=['ros2', 'run', 'plotjuggler', 'plotjuggler', '-l', plotjuggler_config],
            output='screen',
        ),


    ])
