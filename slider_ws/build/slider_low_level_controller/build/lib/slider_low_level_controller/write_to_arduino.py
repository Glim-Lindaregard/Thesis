import rclpy
import serial
from rclpy.node import Node
from std_msgs.msg import UInt8MultiArray

class WriteToArduino(Node):
    def __init__(self):
        super().__init__('write_to_arduino')
        self.declare_parameter('port_name', '/dev/ttyUSB0')
        self.declare_parameter('baud_rate', 9600)

        port = self.get_parameter('port_name').value
        baud = self.get_parameter('baud_rate').value

        self.board = serial.Serial(port, baud, timeout=1)
        self.create_subscription(UInt8MultiArray, '/eight_thrust_pulse', self.callback, 10)

    def callback(self, msg):
        try:
            self.board.write(bytes(msg.data))
        except serial.SerialException as e:
            self.get_logger().error(f'Serial write failed: {e}')

def main(args=None):
    rclpy.init(args=args)
    node = WriteToArduino()
    rclpy.spin(node)
    node.board.close()
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
