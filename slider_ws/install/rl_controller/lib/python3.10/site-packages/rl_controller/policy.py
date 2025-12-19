import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Vector3, PointStamped
from rclpy.qos import QoSPresetProfiles
from nav_msgs.msg import Odometry
from tf_transformations import euler_from_quaternion
from rclpy.qos import QoSProfile, ReliabilityPolicy, DurabilityPolicy, HistoryPolicy
import pprint
from std_msgs.msg import Float64MultiArray

import numpy as np
from numpy.typing import NDArray
import torch 
from stable_baselines3 import PPO
from stable_baselines3.ppo.policies import MultiInputPolicy
from gymnasium import spaces

from time import sleep
import os

MAX_FORCE = 0.7
MAX_TAU = 0.0098

import argparse, rclpy


import numpy as np
from slider_isaac_lab.utils.squashed_ac_sb3 import CustomActorCriticPolicy, DictObsExtractor

class Policy(Node):
    def __init__(self):
        super().__init__('policy')
        
        
        self.declare_parameter('policy_file_path', '')
        self.declare_parameter('frame_id', 'odom')     
        self.declare_parameter('publish_rate', 10.0)  

        self.policy_file_path = self.get_parameter('policy_file_path').value
        self.target_frame = self.get_parameter('frame_id').value
        self.publish_rate = float(self.get_parameter('publish_rate').value)


        self.last_time = self.get_clock().now()
        self.timer = self.create_timer(0.1, self.timer_callback)
        
        self.create_subscription(
            Odometry, 'odom', self.callback, 
            QoSPresetProfiles.get_from_short_key('system_default')
            )
                
        self.target_pt: Odometry | None = None
        self.target_frame = 'odom'  
        
        self.create_subscription(
            Odometry, 'target_point', self.target_callback,
            QoSPresetProfiles.get_from_short_key('system_default')
            )

        self.publisher_ = self.create_publisher(
            Vector3, '/thrust_cmd', 
            QoSPresetProfiles.get_from_short_key('system_default')
            )
        
        self.limits: np.ndarray | None = None  
        
        self.create_subscription(
            Float64MultiArray,
            'action_bounds',
            self.limits_callback,
            QoSPresetProfiles.get_from_short_key('system_default')
        )

        observation_space = spaces.Dict({
            "observations": spaces.Box(low=-np.inf, high=+np.inf, shape=(8,), dtype=np.float32),
            "action_bounds": spaces.Box(low=-np.inf, high=+np.inf, shape=(2,), dtype=np.float32)
        })
        
        action_space = spaces.Box(low=-1, high=1, shape=(2,), dtype=np.float32)
        
        
        policy_kwargs = dict(
        # log_std_init=-1.5,
        features_extractor_class=DictObsExtractor,   # also a class, not string
        features_extractor_kwargs=dict(features_dim=64),
    )

        self.policy = CustomActorCriticPolicy(
            observation_space=observation_space,
            action_space=action_space,
            lr_schedule=lambda progress: 0.001 * progress,
            **policy_kwargs)

        state_dict = torch.load(self.policy_file_path, map_location="cpu")
        self.policy.load_state_dict(state_dict)
        self.policy.eval()  # important for inference

            
            # self.policy = PPO.load(policy_file_path_arg, 
            #                       device="cpu")
        

    def target_callback(self, msg: PointStamped):
        if msg.header.frame_id and msg.header.frame_id != self.target_frame:
            self.get_logger().warn(
                f"target_point in '{msg.header.frame_id}', expected '{self.target_frame}'. "
                "Assuming same frame; transform if needed."
            )
        self.target_pt = msg
        
        
    def timer_callback(self):
        now = self.get_clock().now()
        dt = (now - self.last_time).nanoseconds * 1e-9
        self.last_time = now
        # self.get_logger().info(f"dt_ctrl: {dt:.3f} s")

    def limits_callback(self, msg: Float64MultiArray):
        self.limits = np.array(msg.data, dtype=np.float64).reshape(2, 2)


    def callback(
        self, 
        msg: Odometry
        ):  
        
        if (self.target_pt is None) or (self.limits is None):
            return

        x_ref = float(self.target_pt.pose.pose.position.x)
        y_ref = float(self.target_pt.pose.pose.position.y)
        vx_ref = float(self.target_pt.twist.twist.linear.x)
        vy_ref = float(self.target_pt.twist.twist.linear.y)

        x = float(msg.pose.pose.position.x)
        y = float(msg.pose.pose.position.y)
        vx = float(msg.twist.twist.linear.x)
        vy = float(msg.twist.twist.linear.y)
        omegaz = float(msg.twist.twist.angular.z)
        
        q = msg.pose.pose.orientation
        _, _, theta = euler_from_quaternion([q.x, q.y, q.z, q.w])

        obs = np.array([x,
                        y,
                        vx, 
                        vy,
                        x-x_ref,
                        y-y_ref,
                        vx-vx_ref,
                        vy-vy_ref], dtype=np.float64)
        
                        
        limits = self.limits

        
        obs = {
            "observations":obs,
            "action_bounds":limits
            }
        
        # self.get_logger().info("err=\n%s" % pprint.pformat((x-x_ref, y-y_ref), indent=2, width=80))
        # self.get_logger().info("l2_x=\n%s" % pprint.pformat((x,y), indent=2, width=80))

        forces, _ = self.policy.predict(
            obs, 
            deterministic=True
        )

        f_b = self.world_2_body(fx_w=forces[...,0],
                                fy_w=forces[...,1],
                                theta=theta)
        
        thrust_msg = Vector3()
        thrust_msg.x = float(f_b[0])
        thrust_msg.y = float(f_b[1])
        thrust_msg.z = float(0)
        
        
        
        # self.get_logger().info("forces_raw=\n%s" % pprint.pformat((thrust_msg), indent=2, width=80))

        self.publisher_.publish(thrust_msg)
        
    
    def world_2_body(
        self,
        theta:float,
        fx_w:float,
        fy_w:float
        )->NDArray:
        
        
        c, s = np.cos(theta), np.sin(theta)
        fx_b =  c * fx_w + s * fy_w
        fy_b = -s * fx_w + c * fy_w

        return np.array([fx_b,fy_b])
            
def main(args=None):
    rclpy.init(args=args)
    rl_controller = Policy()
    rclpy.spin(rl_controller)
    rl_controller.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
