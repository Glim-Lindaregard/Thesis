import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Vector3
from rclpy.qos import QoSPresetProfiles
from nav_msgs.msg import Odometry
from tf_transformations import euler_from_quaternion

import numpy as np
import math 

from tsl_optimizer.parameters import MAX_FORCE, MAX_TAU
from slider_thruster_controller.pid_parameters import *





class PIDController(Node):
    def __init__(self):
        super().__init__('pid_controller')
        
        self.declare_parameter('verbose', 1)
        self.declare_parameter('target_x', 1.0)
        self.declare_parameter('target_y', 1.0)
        self.declare_parameter('target_tau', 0.0)
        self.declare_parameter('frequency', 10)

        self.declare_parameter('p', 0.5)
        self.declare_parameter('i', 0.0)
        self.declare_parameter('d', 8.0)

        self.declare_parameter('point_to_point', False)

        self.verbose = self.get_parameter('verbose').get_parameter_value().integer_value
        self.frequency = self.get_parameter('frequency').get_parameter_value().integer_value
        
        self.target_x = self.get_parameter('target_x').get_parameter_value().double_value
        self.target_y = self.get_parameter('target_y').get_parameter_value().double_value
        self.target_tau = self.get_parameter('target_tau').get_parameter_value().double_value
        
        self.p = self.get_parameter('p').get_parameter_value().double_value
        self.i = self.get_parameter('i').get_parameter_value().double_value
        self.d = self.get_parameter('d').get_parameter_value().double_value

        self.point_to_point = self.get_parameter('point_to_point').get_parameter_value().bool_value
        
        self.pub = self.create_publisher(Vector3, '/thrust_cmd', QoSPresetProfiles.get_from_short_key('system_default'))
        self.create_subscription(Odometry, 'odom', self.callback, QoSPresetProfiles.get_from_short_key('system_default'))
        self.create_timer(1/(self.frequency), self.control)

        self.last_sample = self.get_clock().now()
        self.Y = np.array([[0.0], [0.0], [0.0]], dtype=np.float64)
        self.Y_last = self.Y

        self.point_to_point = True 

        if self.point_to_point : 
            self.Yd = np.array([[self.target_x], [self.target_y], [self.target_tau]], dtype=np.float64)
            self.Yd_dot = np.array([[0.0], [0.0], [0.0]], dtype=np.float64)
            self.Yd_dot_dot = np.array([[0.0], [0.0], [0.0]], dtype=np.float64)
            
        
        
    def callback(self, msg: Odometry):
        x = msg.pose.pose.position.x
        y = msg.pose.pose.position.y
        x_q = msg.pose.pose.orientation.x
        y_q = msg.pose.pose.orientation.y
        z_q = msg.pose.pose.orientation.z
        w_q = msg.pose.pose.orientation.w

        euler = euler_from_quaternion([x_q, y_q, z_q, w_q]) # you get roll pitch yaw
        theta = euler[2]

        self.Y_last = self.Y
        self.Y = np.array([[x], [y], [theta]], dtype=np.float64)

        
        
        
    def control(self):
        self.p = self.get_parameter('p').get_parameter_value().double_value
        self.i = self.get_parameter('i').get_parameter_value().double_value
        self.d = self.get_parameter('d').get_parameter_value().double_value

        self.verbose = self.get_parameter('verbose').get_parameter_value().integer_value
        self.frequency = self.get_parameter('frequency').get_parameter_value().integer_value

        self.target_x = self.get_parameter('target_x').get_parameter_value().double_value
        self.target_y = self.get_parameter('target_y').get_parameter_value().double_value
        self.target_tau = self.get_parameter('target_tau').get_parameter_value().double_value

        self.point_to_point = self.get_parameter('point_to_point').get_parameter_value().bool_value
        
        current_time = self.get_clock().now()
        
        if self.point_to_point: 
            self.Yd = np.array([[self.target_x], [self.target_y], [self.target_tau]], dtype=np.float64)
        else : 
            a = 0.25 
            omega = 0.1
            t = current_time.nanoseconds / 1e9
            self.Yd = np.array([[a*math.cos(omega*t)], [a*math.sin(omega*t)], [0.0]], dtype=np.float64)
            self.Yd_dot = np.array([[-a*omega*math.sin(omega*t)], [a*omega*math.cos(omega*t)], [0.0]], dtype=np.float64)
            self.Yd_dot_dot = np.array([[-a*omega*omega*math.cos(omega*t)], [-a*omega*omega*math.sin(omega*t)], [0.0]], dtype=np.float64)

        
        duration = (current_time - self.last_sample).nanoseconds / 1e9

        self.Y_dot = (1/duration) * (self.Y - self.Y_last)

        U = self.PID(self.Y[2], self.Yd_dot_dot, self.Y_dot, self.Yd_dot, self.Yd, self.Y, duration)

        if self.verbose : 
            self.get_logger().info(f'Fx: {U[0]}, Fy: {U[1]}, tau : {U[2]}')

        msg = Vector3()

        msg.x = float(U[0])
        msg.y = float(U[1])
        msg.z = float(U[2])


        self.last_sample = current_time

        self.pub.publish(msg)
        

    def PID(self, theta, Yd_dot_dot, Y_dot, Yd_dot, Yd, Y, duration): 
        self.G = np.array(
            [[(1/mass) * math.cos(theta),     -(1/mass) * math.sin(theta),      0],
             [(1/mass) * math.sin(theta),     (1/mass) * math.cos(theta),       0],
             [0,                            0,                          1/Izz]], 
            dtype=np.float64)

        G_inv = np.linalg.inv(self.G)

        e_dot = Y_dot - Yd_dot 
        e = Y - Yd 

        U_new = G_inv.dot(Yd_dot_dot - self.d * gain_K1.dot(e_dot) - self.p * gain_K2.dot(e) - self.i * duration * e )

        return self.constrain(U_new)
    
    def constrain(self, U): 
        fx = U[0]
        fy = U[1]
        tau = U[2]

        fx = min(max(fx, -MAX_FORCE), MAX_FORCE)
        fy = min(max(fy, -MAX_FORCE), MAX_FORCE)
        tau = min(max(tau, -MAX_TAU), MAX_TAU)

        U[0] = fx 
        U[1] = fy 
        U[2] = tau 

        return U 

        
        
    
def main(args=None):
    rclpy.init(args=args)
    node = PIDController()
    rclpy.spin(node)        
    node.destroy_node()
    rclpy.shutdown()
            

if __name__ == '__main__':
    main()