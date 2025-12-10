import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Vector3
from std_msgs.msg import UInt8MultiArray
from rclpy.qos import QoSPresetProfiles

import numpy as np

from tsl_optimizer.solver import Solver
from tsl_optimizer.parameters import MAX_FORCE



class ThrustController(Node):
    def __init__(self):
        super().__init__('pwm_publisher')
        
        self.declare_parameter('verbose_optimizer', 0)
        self.declare_parameter('verbose_pwm', 0)
        self.declare_parameter('frequency', 10)
        self.declare_parameter('resolution', 64)

        self.verbose_optimizer = self.get_parameter('verbose_optimizer').get_parameter_value().integer_value
        self.verbose_pwm = self.get_parameter('verbose_pwm').get_parameter_value().integer_value
        self.frequency = self.get_parameter('frequency').get_parameter_value().integer_value
        self.resolution = self.get_parameter('resolution').get_parameter_value().integer_value

        self.solver = Solver()
        
        self.signals = [self.create_pwm(0, self.resolution) for i in range(8)]
        self.i = 0
        
        self.create_subscription(Vector3, 'thrust_cmd', self.callback, QoSPresetProfiles.get_from_short_key('system_default'))
        self.pub = self.create_publisher(UInt8MultiArray, '/slider_controller', QoSPresetProfiles.get_from_short_key('system_default'))

        self.create_timer(1/(self.frequency * self.resolution), self.send_signals)
        
        
    def callback(self, msg: Vector3):
        U = np.array([[msg.x], [msg.y], [msg.z]], dtype=np.float64)
        T = self.solver.solve(U)
        if self.verbose_optimizer > 0: 
            self.get_logger().info(f'\n Fx = {msg.x: 2.2f}\n Fy = {msg.y: 2.2f}\ntau = {msg.z: 2.2f}')
            self.get_logger().info(f'cmd: {T}')

        self.signals = [self.create_pwm(T[i], self.resolution) for i in range(8)]
        for n, s in enumerate(self.signals):
            if s[0] == 1 and s[1] == 0:
                self.signals[n][1] = 1
        
        
    def send_signals(self):
        self.verbose_optimizer = self.get_parameter('verbose_optimizer').get_parameter_value().integer_value
        self.verbose_pwm = self.get_parameter('verbose_pwm').get_parameter_value().integer_value
        self.frequency = self.get_parameter('frequency').get_parameter_value().integer_value

        req = UInt8MultiArray()
        req.data = [int(self.signals[i][self.i]) for i in range(0, 8)]

        self.i += 1
        self.i %= self.resolution
        
        self.pub.publish(req)

    def create_pwm(self, thrust, resolution):
        if self.verbose_pwm > 0:    
            self.get_logger().info(f"thrust value to be modulated : {thrust}")
        value = thrust / MAX_FORCE
        unit = MAX_FORCE / resolution
        
        if value < unit:
            value = 0.0
        
        
        number_of_pulses = int(value /unit)

        signals = [1 for _ in  range(number_of_pulses)]
        signals+= [0 for _ in  range(resolution - number_of_pulses)] 
        if self.verbose_pwm > 0:
            self.get_logger().info(f"pwm signals : {signals}")
        return signals
        
        
    
def main(args=None):
    rclpy.init(args=args)
    node = ThrustController()
    rclpy.spin(node)        
    node.destroy_node()
    rclpy.shutdown()
            

if __name__ == '__main__':
    main()