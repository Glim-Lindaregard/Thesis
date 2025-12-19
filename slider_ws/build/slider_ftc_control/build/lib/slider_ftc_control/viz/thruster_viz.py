#!/usr/bin/env python3
import rclpy
from rclpy.node import Node
from std_msgs.msg import UInt8MultiArray

import matplotlib.pyplot as plt
import numpy as np

print("Starting thruster visualization node...")
class ThrusterViz(Node):
    def __init__(self):
        super().__init__('thruster_viz')

        self.sub = self.create_subscription(
            UInt8MultiArray,
            'eight_thrust_pulse',
            self.callback,
            10
        )

        # latest pulse bits (0/1) for 8 thrusters
        self.pulses = [0] * 8

        # set up figure
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal')
        self.ax.set_xlim(-1.5, 1.5)
        self.ax.set_ylim(-1.5, 1.5)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title('Thruster Pulse Visualization')

        # draw slider body (square)
        body = plt.Rectangle((-1, -1), 2, 2, fill=False, linewidth=2)
        self.ax.add_patch(body)

        # thruster positions (around corners)
        # order: 0..7
        self.positions = np.array([
            [ 0.9, -0.9],  # T1 (bottom-right)
            [ 0.9, -1.1],  # T2
            [-0.9, -0.9],  # T3 (bottom-left)
            [-0.9, -1.1],  # T4
            [-0.9,  0.9],  # T5 (top-left)
            [-0.9,  1.1],  # T6
            [ 0.9,  0.9],  # T7 (top-right)
            [ 0.9,  1.1],  # T8
        ])

        # simple arrow directions (you can tweak to match real layout)
        self.directions = np.array([
            [ 0.3,  0.0],  # T1
            [ 0.0,  0.3],  # T2
            [-0.3,  0.0],  # T3
            [ 0.0,  0.3],  # T4
            [-0.3,  0.0],  # T5
            [ 0.0, -0.3],  # T6
            [ 0.3,  0.0],  # T7
            [ 0.0, -0.3],  # T8
        ])

        # create arrow artists
        self.arrows = []
        for (x, y), (dx, dy) in zip(self.positions, self.directions):
            arr = self.ax.arrow(
                x, y, dx, dy,
                head_width=0.08,
                head_length=0.08,
                length_includes_head=True,
                color='gray'
            )
            self.arrows.append(arr)

        # enable interactive mode
        plt.ion()
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def callback(self, msg: UInt8MultiArray):
        # store latest pulses (clamp len to 8)
        data = list(msg.data)
        if len(data) >= 8:
            self.pulses = data[:8]
        else:
            # pad if shorter
            self.pulses = data + [0] * (8 - len(data))

    def update_plot(self):
        # update arrow colors based on pulses
        for i, arr in enumerate(self.arrows):
            active = bool(self.pulses[i])
            color = 'green' if active else 'gray'
            arr.set_color(color)

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()


def main():
    rclpy.init()
    node = ThrusterViz()

    try:
        # manual spin/update loop so matplotlib and ROS can coexist
        rate_hz = 50.0
        dt = 1.0 / rate_hz
        while rclpy.ok():
            rclpy.spin_once(node, timeout_sec=0.0)
            node.update_plot()
            plt.pause(dt)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()
        plt.close('all')


if __name__ == '__main__':
    main()
