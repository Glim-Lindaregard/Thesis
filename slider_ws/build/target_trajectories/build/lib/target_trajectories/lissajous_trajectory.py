#!/usr/bin/env python3
import math
from typing import Tuple
import rclpy
from rclpy.node import Node
from nav_msgs.msg import Odometry
from rclpy.qos import QoSPresetProfiles

def lissajous_state(
    t: float,
    center: Tuple[float, float, float],
    A: float, B: float,
    a: int, b: int,
    delta: float,
    omega: float,
):
    cx, cy, cz = center
    theta = omega * t
    x = cx + A * math.sin(a * theta + delta)
    y = cy + B * math.sin(b * theta)
    z = cz
    vx = A * a * omega * math.cos(a * theta + delta)
    vy = B * b * omega * math.cos(b * theta)
    vz = 0.0
    return (x, y, z), (vx, vy, vz)

class LissajousTrajectory(Node):
    def __init__(self):
        super().__init__('lissajous_trajectory')

        self.declare_parameter('frame_id', 'world')
        self.declare_parameter('publish_rate', 10.0)
        self.declare_parameter('lissa_A', 2.0)
        self.declare_parameter('lissa_B', 2.0)
        self.declare_parameter('lissa_a', 1)
        self.declare_parameter('lissa_b', 2)
        self.declare_parameter('lissa_delta', 0.0)
        self.declare_parameter('lissa_omega', 0.05)

        frame = self.get_parameter('frame_id').value
        rate = float(self.get_parameter('publish_rate').value)
        self.A = float(self.get_parameter('lissa_A').value)
        self.B = float(self.get_parameter('lissa_B').value)
        self.a = int(self.get_parameter('lissa_a').value)
        self.b = int(self.get_parameter('lissa_b').value)
        self.delta = float(self.get_parameter('lissa_delta').value)
        self.omega = float(self.get_parameter('lissa_omega').value)

        self.dt = 1.0 / rate
        self.t = 0.0

        self.msg = Odometry()
        self.msg.header.frame_id = frame
        (x, y, z), (vx, vy, vz) = lissajous_state(
            0.0, (0.0, 0.0, 0.0),
            self.A, self.B, self.a, self.b, self.delta, self.omega,
        )
        self.msg.pose.pose.position.x = x
        self.msg.pose.pose.position.y = y
        self.msg.pose.pose.position.z = z
        self.msg.twist.twist.linear.x = vx
        self.msg.twist.twist.linear.y = vy
        self.msg.twist.twist.linear.z = vz

        self.publisher_ = self.create_publisher(
            Odometry, 'target_point',
            QoSPresetProfiles.get_from_short_key('system_default')
        )
        self.timer = self.create_timer(self.dt, self._tick)

    def _tick(self):
        self.t += self.dt
        (x, y, z), (vx, vy, vz) = lissajous_state(
            self.t, (0.0, 0.0, 0.0),
            self.A, self.B, self.a, self.b, self.delta, self.omega,
        )
        self.msg.header.stamp = self.get_clock().now().to_msg()
        self.msg.pose.pose.position.x = x
        self.msg.pose.pose.position.y = y
        self.msg.pose.pose.position.z = z
        self.msg.twist.twist.linear.x = vx
        self.msg.twist.twist.linear.y = vy
        self.msg.twist.twist.linear.z = vz
        self.publisher_.publish(self.msg)

def main(args=None):
    rclpy.init(args=args)
    node = LissajousTrajectory()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
