import rclpy
import numpy as np
from rclpy.node import Node
from nav_msgs.msg import Odometry
from tf_transformations import euler_from_quaternion
from std_msgs.msg import UInt8MultiArray
from rclpy.qos import QoSPresetProfiles
from slider_ftc_control.config import make_default_config
from slider_ftc_control.ams import build_ams_cache, allocate_wrench
from slider_ftc_control.ams.cache import mask_to_index
from slider_ftc_control.mpc import MPCController,compute_wrench_bounds
from tf_transformations import quaternion_matrix




# ------------------- SET FAILED THRUSTERS HERE (Outside FTCNode so that the plotting function can reach it) --------------------
START_HEALTH = [1,1,1,1,1,1,1,1]  # [T11, T12, T21, T22, T31, T32, T41, T42]   You can put 0 here to fail a thruster from t=0
FAILURE_TIMES  = [None,  None,  None, None, None,    None,    None, None]   # What time to inject failure (None to never fail)
FAILURE_STATES = [2,     2,    2,    2,    1,    2,    2,    2] # 0 -> Passive (turned off), 1-> healthy, 2-> active (always on)
# -------------------------------------------------------------------------------------------------------------------------------


class FTCNode(Node):
    def __init__(self):
        super().__init__('ftc_node')

        self.health = START_HEALTH.copy()
        self.failure_times  = FAILURE_TIMES.copy()
        self.failure_states = FAILURE_STATES.copy()

        self.cfg = make_default_config()
        self.ams_cache = build_ams_cache(self.cfg)
        self.mpc = MPCController(self.cfg)

        self.state = None
        self.reference = None

        # --- Stuff for PWM ---
        self.declare_parameter('pwm_frequency', 10)
        self.declare_parameter('pwm_resolution', 64)
        self.declare_parameter('max_force', 0.7)

        self.pwm_frequency = self.get_parameter('pwm_frequency').get_parameter_value().integer_value
        self.pwm_resolution = self.get_parameter('pwm_resolution').get_parameter_value().integer_value
        self.max_force = self.get_parameter('max_force').get_parameter_value().double_value
        
        self.signals = [self.create_pwm(0.0) for _ in range(self.cfg.phys.N_thrusters)]
        self.i = 0

        # -PWM output timer -
        self.pwm_timer = self.create_timer(
            1.0 / (self.pwm_frequency * self.pwm_resolution),
            self.send_signals,
        )


        # --- control refresh time ---
        self.refresh_rate = 0.1
        self.timer = self.create_timer(self.refresh_rate, self.control_step)


        # --- This is just to slow down the printing of states and stuff in the terminal ---
        self.state_print_counter = 0
        self.ref_print_counter = 0
        self.ad_print_counter = 0
        self.printer_timing = 3


        self.failure_already_logged = [False]*self.cfg.phys.N_thrusters

        
        # --- Define subscribers and publishers ---

        # - \odom sub -
        self.state_sub = self.create_subscription(
            Odometry,
            '/odom',
            self.odom_callback,
            QoSPresetProfiles.get_from_short_key('system_default')
        )

        # - \target_point sub - 
        self.ref_sub = self.create_subscription(
            Odometry,
            '/target_point',
            self.target_point_callback,
            QoSPresetProfiles.get_from_short_key('system_default')
        )

        # - /eight_thrust_pulse pub -
        self.cmd_pub = self.create_publisher(
            UInt8MultiArray,
            '/eight_thrust_pulse',
            QoSPresetProfiles.get_from_short_key('system_default')
        )


        self.get_logger().info('FTCNode started, waiting for /odom and /target_point...')

        self.t = 0.0 # For failure activation 

        
    def odom_callback(self, msg: Odometry):

        x = float(msg.pose.pose.position.x)
        y = float(msg.pose.pose.position.y)


        

        q_msg = msg.pose.pose.orientation
        q = np.array([q_msg.x, q_msg.y, q_msg.z, q_msg.w], dtype=float)
        _, _, theta = euler_from_quaternion([q[0], q[1], q[2], q[3]])
        vx_b = float(msg.twist.twist.linear.x)
        vy_b = float(msg.twist.twist.linear.y)

        # Rotation from body -> world (from the pose quaternion)
        R = quaternion_matrix(q)[:3, :3]

        v_world = R @ np.array([vx_b, vy_b, 0.0]).reshape(3)
        vx = float(v_world[0])
        vy = float(v_world[1])
        
        omega = float(msg.twist.twist.angular.z)


        self.state = np.array([x, y, q[0], q[1], q[2], q[3], vx, vy, omega], dtype=float)

        self.state_print_counter += 1
        if self.state_print_counter % self.printer_timing == 0:
            self.get_logger().info(
                f"STATE: x={x:.2f}, y={y:.2f}, th={theta:.2f}, "
                f"vx={vx:.2f}, vy={vy:.2f}, w={omega:.2f} \n"
            )



    def target_point_callback(self, msg: Odometry):

        x_ref = float(msg.pose.pose.position.x)
        y_ref = float(msg.pose.pose.position.y)

        q_msg = msg.pose.pose.orientation
        q = np.array([q_msg.x, q_msg.y, q_msg.z, q_msg.w], dtype=float)
        _, _, theta_ref = euler_from_quaternion([q[0], q[1], q[2], q[3]])


        vx_ref = float(msg.twist.twist.linear.x)
        vy_ref = float(msg.twist.twist.linear.y)
        omega_ref = float(msg.twist.twist.angular.z)
        

        self.reference = np.array(
            [x_ref, y_ref, q[0], q[1], q[2], q[3], vx_ref, vy_ref, omega_ref],
            dtype=float
        )

        self.ref_print_counter += 1
        if self.ref_print_counter % self.printer_timing == 0:
            self.get_logger().info(
                f"REF: x={x_ref:.2f}, y={y_ref:.2f}, th={theta_ref:.2f}, "
                f"vx={vx_ref:.2f}, vy={vy_ref:.2f}, w={omega_ref:.2f} \n"
            )


    def control_step(self):
        if self.state is None:
            self.get_logger().warning('State or reference not yet available for control step.')
            return

        self.t += self.refresh_rate

        health_arr = np.array(self.health, dtype=int)

        # --- Print message if thruster failes ---
        RED = "\033[91m"
        BOLD = "\033[1m"
        RESET = "\033[0m"

        for i in range(self.cfg.phys.N_thrusters):
            t_fail = self.failure_times[i]
            if t_fail is not None and self.t >= t_fail:
                health_arr[i] = int(self.failure_states[i])
                if not self.failure_already_logged[i]:
                    msg = (
                        f"{RED}{BOLD}"
                        f"*** THRUSTER {i+1} FAILED (code {self.failure_states[i]}) ***"
                        f"{RESET}"
                    )
                    self.get_logger().error(msg)
                    self.failure_already_logged[i] = True


        #--------------- Check if too many failed thrusters------------------#
        ams_mask = (health_arr != 0).astype(int)
        n_thrusters = self.cfg.phys.N_thrusters
        n_failed = n_thrusters - int(ams_mask.sum())
        if n_failed > self.cfg.phys.max_failed_thr:
            self.get_logger().error(
                f"Too many failed thrusters: {n_failed} > max_failed_thr={self.cfg.phys.max_failed_thr}. "
                "Controller is not designed for this failure level. Shutting down."
            )
            rclpy.shutdown()
            raise SystemExit(1)
        #-------------------------------------------------------------------#

        idx = mask_to_index(health_arr)
        mode = self.ams_cache[idx]

        A_curr = mode.A_mode

        u_min = self.cfg.phys.u_min
        u_max = self.cfg.phys.u_max
        
        # --- MPC bounds for current failure case
        bounds = compute_wrench_bounds(
            A=A_curr,
            u_min=u_min,
            u_max=u_max,
            )

        # --- Main control step ---
        a_d = self.mpc.step(
            x_current = self.state, 
            x_ref = self.reference,
            bounds = bounds
        )

        # --- Allocation ---
        try:
            cmd = allocate_wrench(
                a_d=a_d,
                mode=mode,
                cfg=self.cfg,
            )
        except ValueError as e:
            self.get_logger().warn(f"Allocation failed: {e}. Setting all thrusters to 0.")
            cmd = np.zeros(self.cfg.phys.N_thrusters, dtype=float)


        # --- After step, set passive to 0 and active to u_max
        for i, h in enumerate(health_arr):
            if h == 0:
                cmd[i] = 0.0
            elif h == 2:
                cmd[i] = self.cfg.phys.u_max[i]

        # --- Convert continous 'cmd' to PWM'ed 'signals'
        self.signals = [self.create_pwm(float(cmd[i])) for i in range(self.cfg.phys.N_thrusters)]


        # --- Start spike protection form other PWM code (CHECK WITH NEKTARIOS IF THIS IS NECESSARY) ---
        for n, s in enumerate(self.signals):
            if len(s) >= 2 and s[0] == 1 and s[1] == 0:
                self.signals[n][1] = 1

        # --- Print commands ---
        self.ad_print_counter +=1
        if self.ad_print_counter % self.printer_timing == 0:
           self.get_logger().info(
                  f"a_d = [{a_d[0]:.2f}, {a_d[1]:.2f}, {a_d[2]:.2f}], "
                  f"u = [{', '.join(f'{x:.2f}' for x in cmd)}] \n"
)



    def create_pwm(self, thrust: float) -> list[int]:

        duty = max(0.0, min(1.0, thrust / self.max_force))

        number_of_pulses = int(duty * self.pwm_resolution)

        signals = [1] * number_of_pulses + [0] * (self.pwm_resolution - number_of_pulses)
        return signals


    def send_signals(self):

        msg = UInt8MultiArray()
        msg.data = [int(self.signals[i][self.i]) for i in range(self.cfg.phys.N_thrusters)]

        self.i = (self.i + 1) % self.pwm_resolution

        self.cmd_pub.publish(msg)



def main(args=None):
    rclpy.init(args=args)
    node = FTCNode()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        node.get_logger().info('Shutting down FTCNode (Ctrl-C).')
    finally:
        node.destroy_node()
        if rclpy.ok():
            rclpy.shutdown()


if __name__ == '__main__':
    main()
