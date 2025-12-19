# target_publisher.py
import rclpy
from rclpy.node import Node
from nav_msgs.msg import Odometry     # changed
from rclpy.qos import QoSPresetProfiles
import numpy as pn
from tf_transformations import quaternion_from_euler

class TargetPoint(Node):
    def __init__(self):
        super().__init__('target_point')

        self.declare_parameter('target_x', 0.0)
        self.declare_parameter('target_y', 0.0)
        self.declare_parameter('target_theta',0.0)
    
        self.declare_parameter('frame_id', 'world')
        self.declare_parameter('publish_rate', 10.0)  
                
        frame = self.get_parameter('frame_id').get_parameter_value().string_value
        rate = self.get_parameter('publish_rate').get_parameter_value().double_value
        self.target_x = float(self.get_parameter('target_x').value)
        self.target_y = float(self.get_parameter('target_y').value)
        self.target_theta = float(self.get_parameter('target_theta').value)

        self.msg = Odometry()
        self.msg.header.frame_id = frame
        self.msg.pose.pose.position.x = self.target_x
        self.msg.pose.pose.position.y = self.target_y
    

        qx, qy, qz, qw = quaternion_from_euler(0.0, 0.0, self.target_theta)

        self.msg.pose.pose.orientation.x = qx
        self.msg.pose.pose.orientation.y = qy
        self.msg.pose.pose.orientation.z = qz
        self.msg.pose.pose.orientation.w = qw


        self.msg.twist.twist.linear.x = 0.0
        self.msg.twist.twist.linear.y = 0.0
        self.msg.twist.twist.linear.z = 0.0

        self.publisher_ = self.create_publisher(Odometry, 'target_point', 
                                                QoSPresetProfiles.get_from_short_key('system_default'))
        self.timer = self.create_timer(1.0 / rate, self._tick)

    def _tick(self):
        self.msg.header.stamp = self.get_clock().now().to_msg()
        self.publisher_.publish(self.msg)

def main(args=None):
    rclpy.init(args=args)
    node = TargetPoint()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()
