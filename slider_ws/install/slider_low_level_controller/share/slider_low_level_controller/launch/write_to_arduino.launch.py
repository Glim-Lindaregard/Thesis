from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration

def generate_launch_description():
    return LaunchDescription([
        DeclareLaunchArgument('port_name', default_value='/dev/ttyUSB0'),
        DeclareLaunchArgument('baud_rate', default_value='9600'),
        Node(
            package='slider_low_level_controller',
            executable='write_to_arduino',
            parameters=[{
                'port_name': LaunchConfiguration('port_name'),
                'baud_rate': LaunchConfiguration('baud_rate')
            }]
        )
    ])
