import rclpy
from rclpy.node import Node
from std_msgs.msg import UInt8MultiArray
from visualization_msgs.msg import Marker, MarkerArray
from geometry_msgs.msg import Point

class ThrusterMarkers(Node):
    def __init__(self):
        super().__init__('thruster_markers')

        # Thruster positions in BODY frame [m].
        # Order must match your A-matrix / real thruster order.
        # Example from the Slider paper: positions in mm, convert to m.
        self.thr_pos = [
            ( 0.140,  0.195),  # T1
            ( 0.195,  0.140),  # T2
            ( 0.140, -0.195),  # T3
            ( 0.195, -0.140),  # T4
            (-0.140, -0.195),  # T5
            (-0.195, -0.140),  # T6
            (-0.140,  0.195),  # T7
            (-0.195,  0.140),  # T8
        ]


        # Rough thrust directions in BODY frame (change to match your A-matrix)
        self.thr_dir = [
            ( 0.0,  1.0),   # T1
            ( 1.0,  0.0),   # T2
            ( 0.0, -1.0),   # T3
            ( 1.0,  0.0),   # T4
            ( 0.0, -1.0),   # T5
            (-1.0, 0.0),    # T6
            ( 0.0,  1.0),   # T7
            (-1.0, 0.0),    # T8
        ]



        self.sub = self.create_subscription(
            UInt8MultiArray,
            '/thruster_state',
            self.cb_state,
            10,
        )

        self.pub = self.create_publisher(
            MarkerArray,
            '/thruster_markers',
            10,
        )

        self.get_logger().info('ThrusterMarkers node started.')

    def cb_state(self, msg: UInt8MultiArray):
        states = list(msg.data)
        n = len(self.thr_pos)
        if len(states) < n:
            states += [0] * (n - len(states))

        ma = MarkerArray()
        now = self.get_clock().now().to_msg()

        arrow_len = 0.15  # [m] length of arrow

        for i, ((x, y), (dx, dy)) in enumerate(zip(self.thr_pos, self.thr_dir)):
            m = Marker()
            m.header.frame_id = 'body'   # keep body; it will move with the slider if TF is set up
            m.header.stamp = now
            m.ns = 'thrusters'
            m.id = i
            m.type = Marker.ARROW
            m.action = Marker.ADD

            # Arrow defined by two points: start at thruster, end in thrust direction
            start = Point()
            start.x = float(x)
            start.y = float(y)
            start.z = 0.05

            end = Point()
            end.x = float(x + arrow_len * dx)
            end.y = float(y + arrow_len * dy)
            end.z = 0.05

            m.points = [start, end]

            # For ARROW with points:
            # scale.x = shaft diameter, scale.y = head diameter, scale.z = head length
            m.scale.x = 0.02  # shaft thickness
            m.scale.y = 0.05  # head size
            m.scale.z = 0.05  # head length

            state = states[i]

            # Colors: idle = grey, firing = green, broken = red
            if state == 0:
                m.color.r = 0.5
                m.color.g = 0.5
                m.color.b = 0.5
            elif state == 1:
                m.color.r = 0.0
                m.color.g = 1.0
                m.color.b = 0.0
            else:  # state == 2
                m.color.r = 1.0
                m.color.g = 0.0
                m.color.b = 0.0
            m.color.a = 0.9

            m.lifetime.sec = 0
            ma.markers.append(m)

        self.pub.publish(ma)



def main(args=None):
    rclpy.init(args=args)
    node = ThrusterMarkers()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
