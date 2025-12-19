from launch import LaunchDescription
from launch_ros.actions import Node 
from ament_index_python.packages import get_package_share_directory
from launch.substitutions import LaunchConfiguration

import xacro 
import os 

def generate_launch_description():

    position = [0.0, 0.0, 0.0]
    
    # [Roll, Pitch, Yaw]
    orientation = [0.0, 0.0, 0.0]    

    entity_name = "slider"
    
    xacro_file = "slider.urdf.xacro"

    package_description = "slider_description"

    use_sim_time = LaunchConfiguration('use_sim_time', default='true')

    robot_desc_path = os.path.join(get_package_share_directory(package_description), "urdf", xacro_file)
    urdf = xacro.process_file(robot_desc_path).toxml()

    publish_robot_description = Node(
        package='slider_gazebo',
        executable='urdf_publisher.py',
        name='urdf_publisher',
        output='screen',
        arguments=['-xml_string', urdf,
                   '-robot_description_topic', '/robot_description'
                   ]
    )


    robot_state_publisher_node = Node(
        package='robot_state_publisher',
        executable='robot_state_publisher',
        parameters=[{'use_sim_time': use_sim_time, 'robot_description': urdf}],
        output="screen"
    )
    

    joint_state_publisher_node = Node(
            package='joint_state_publisher',
            executable='joint_state_publisher',
            name='joint_state_publisher',
            output='screen',
    )



    spawn_robot = Node(
        package='gazebo_ros',
        executable='spawn_entity.py',
        name='spawn_entity',
        output='screen',
        arguments=['-entity',
                   entity_name,
                   '-x', str(position[0]), '-y', str(position[1]
                                                     ), '-z', str(position[2]),
                   '-R', str(orientation[0]), '-P', str(orientation[1]
                                                        ), '-Y', str(orientation[2]),
                   '-topic', '/robot_description'
                   ]
    )


    return LaunchDescription([
        publish_robot_description, 
        joint_state_publisher_node, 
        robot_state_publisher_node,
        spawn_robot
    ])