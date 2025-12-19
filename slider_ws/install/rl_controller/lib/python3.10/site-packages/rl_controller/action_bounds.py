import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Vector3, PointStamped
from rclpy.qos import QoSPresetProfiles
from nav_msgs.msg import Odometry
from tf_transformations import euler_from_quaternion
from rclpy.qos import QoSProfile, ReliabilityPolicy, DurabilityPolicy, HistoryPolicy
import pprint

import numpy as np
from numpy.typing import NDArray
import torch 
from stable_baselines3 import PPO
from stable_baselines3.ppo.policies import MultiInputPolicy
from gymnasium import spaces
# from slider_gym.gyms.go_2_point_gym import get_allowable_control
from std_msgs.msg import Float64MultiArray, MultiArrayDimension, MultiArrayLayout


from slider_isaac_lab.utils.get_admissible_control import get_admissible_action_bounds


class ActionBounds(Node):
    def __init__(self):
        super().__init__('action_bounds')
        
        self.create_subscription(
            Odometry, 'odom', self.callback, 
            QoSPresetProfiles.get_from_short_key('system_default')
            )
            
        self.limits_pub = self.create_publisher(
            Float64MultiArray, 'action_bounds',
            QoSPresetProfiles.get_from_short_key('system_default')
            )

        self.declare_parameter('x_bounds', [float(-2), float(2)])
        self.declare_parameter('y_bounds', [float(-2), float(2)])

    
    def callback(self, msg: Odometry):  

        x = float(msg.pose.pose.position.x)
        y = float(msg.pose.pose.position.y)
        vx = float(msg.twist.twist.linear.x)
        vy = float(msg.twist.twist.linear.y)
        omegaz = float(msg.twist.twist.angular.z)
        
        q = msg.pose.pose.orientation
        _, _, theta = euler_from_quaternion([q.x, q.y, q.z, q.w])

        x_bounds_param = self.get_parameter('x_bounds').get_parameter_value().double_array_value
        y_bounds_param = self.get_parameter('y_bounds').get_parameter_value().double_array_value
        x_bounds = np.array(list(x_bounds_param), dtype=np.float64)
        y_bounds = np.array(list(y_bounds_param), dtype=np.float64)

        limits = get_admissible_action_bounds(
            x=torch.Tensor([x, y, theta, vx, vy, omegaz]),
            x_bounds=torch.Tensor([x_bounds, y_bounds])
        )
        
        # limits = np.array(([-0.7,+0.7],
        #                   [-0.7,+0.7]))
        
        limits = np.array(limits)

        limits_mat = np.asarray(limits, dtype=np.float64).reshape(2, 2)
        msg_out = Float64MultiArray()
        msg_out.data = limits_mat.ravel(order='C').tolist()  
        self.limits_pub.publish(msg_out)


def main(args=None):
    rclpy.init(args=args)
    node = ActionBounds()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()
