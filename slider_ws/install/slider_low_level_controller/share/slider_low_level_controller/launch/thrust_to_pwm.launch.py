# slider_thruster.launch.py
#
# ROS 2 Python launch file converted from your XML.
# It declares the same arguments and passes them as typed parameters.

from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node
from launch_ros.parameter_descriptions import ParameterValue

def generate_launch_description():
    
    frequency_arg = DeclareLaunchArgument(
        'frequency', 
        default_value='10',
        description='PWM frequency (Hz)'
    )
    resolution_arg = DeclareLaunchArgument(
        'resolution', 
        default_value='64',
        description='PWM resolution'
    )
    
    max_force_arg = DeclareLaunchArgument(
        'max_force', 
        default_value='0.7',
        description='Max thruster force'
    )


    # Substitute values from command line or defaults
    frequency = LaunchConfiguration('frequency')
    resolution = LaunchConfiguration('resolution')
    max_force = LaunchConfiguration('max_force')


    thrust_to_pwm = Node(
        package='slider_low_level_controller',
        executable='thrust_to_pwm',
        name='thrust_to_pwm',
        output='screen',
        parameters=[{
            'frequency': frequency,
            'resolution': resolution,
            'max_force': max_force,
        }]
    )



    return LaunchDescription([
        frequency_arg,
        resolution_arg,
        max_force_arg,
        thrust_to_pwm,
    ])
