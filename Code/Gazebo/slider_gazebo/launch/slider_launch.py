import os

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import IncludeLaunchDescription
from launch.launch_description_sources import PythonLaunchDescriptionSource


def generate_launch_description():
    pkg_slider_gazebo = get_package_share_directory('slider_gazebo')


    start_world = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(pkg_slider_gazebo, 'launch',
                         'world_launch.py'),
        )
    )


    spawn_robot_world = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(pkg_slider_gazebo, 'launch',
                         'spawn_robot_launch.py'),
        )
    )

    # start_rviz = IncludeLaunchDescription(
    #     PythonLaunchDescriptionSource(
    #         os.path.join(pkg_slider_gazebo, 'launch', 'rviz_launch.py'),
    #     )
    # )

    return LaunchDescription([
        start_world,
        spawn_robot_world
    ])
