import rclpy
from rclpy.node import Node
from nav_msgs.msg import Odometry


class FTCNode(Node):
    def __init__(self):
        super().__init__('ftc_node')

        # Subscribe to /odom
        self.subscription = self.create_subscription(
            Odometry,
            '/odom',
            self.odom_callback,
            10
        )

        self.get_logger().info('FTCNode started, waiting for /odom...')

    def odom_callback(self, msg: Odometry):
        # Extract pose just to prove it works
        x = msg.pose.pose.position.x
        y = msg.pose.pose.position.y
        self.get_logger().info(f'Received odom: x = {x:.3f}, y = {y:.3f}')


def main(args=None):
    rclpy.init(args=args)
    node = FTCNode()
    try:
        rclpy.spin(node)
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
