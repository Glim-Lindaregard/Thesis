# minimal_plot_odom_target_pos_vel.py
# - Reads /odom (x,y,vx,vy) and optional /target_point from a rosbag2
# - Reads /action_bounds (3x2 or 2x2; uses first 2x2 for x/y min/max)
# - Draws ellipses centered on odom points; axes = (xmax-xmin, ymax-ymin) * SCALE
#   (no center shift), and only if bounds != ±DEFAULT_PM on BOTH axes
# - Plot 1: XY trajectory with time colormap, target TRAJECTORY (continuous) + latest target marker, and square [-2,2]
# - Plot 2: Position vs time (x(t), y(t)) with ±2 guide lines
# - Plot 3: Velocity vs time (vx(t), vy(t)) with ±2 guide lines
# - Plot 4: Thrust commands fx, fy (ZOH) with per-sample action-bounds min/max overlays (2-row subplot)
# - Plot 5: Target vs Odom and Error (3-row subplot: x, y, L2), target sampled ZOH at odom timestamps

import argparse
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message
import rosbag2_py

DEFAULT_BAG = "/root/slider_ws/test1/test1_0.db3"
odom_topic = "/odom"
target_topic = "/target_point"         
action_bounds_topic = "/action_bounds" # Float32MultiArray (3x2 or 2x2)
thrust_cmd_topic = "/thrust_cmd"       # geometry_msgs/Vector3: [fx, fy, tau]
DEFAULT_STORAGE_ID = "sqlite3"        

SCALE = 0          
DEFAULT_PM = 0.7   
TOL = 1e-3         

parser = argparse.ArgumentParser(description="Plot odometry/target/commands from a rosbag2.")
parser.add_argument("bag", 
                    nargs="?",
                    default=DEFAULT_BAG,
                    help="Path to rosbag2 file (.db3/.mcap). Defaults to example path.")
parser.add_argument("--storage-id", 
                    default=DEFAULT_STORAGE_ID, 
                    choices=["sqlite3", "mcap"],
                    help="rosbag2 storage id (default: sqlite3).")

args = parser.parse_args()
bag = args.bag
storage_id = args.storage_id

reader = rosbag2_py.SequentialReader()
reader.open(
    rosbag2_py.StorageOptions(uri=bag, storage_id=storage_id),
    rosbag2_py.ConverterOptions("cdr", "cdr")
)

types = {t.name: t.type for t in reader.get_all_topics_and_types()}
if odom_topic not in types:
    raise SystemExit(f"Missing {odom_topic}. Found topics: {list(types)}")

Odom = get_message(types[odom_topic])
Target = get_message(types[target_topic]) if target_topic in types else None
ActionBounds = get_message(types[action_bounds_topic]) if action_bounds_topic in types else None
Thrust = get_message(types[thrust_cmd_topic]) if thrust_cmd_topic in types else None

# Odom series
times, xs, ys, vxs, vys = [], [], [], [], []
ell_bounds = []            

txs, tys = [], []
latest_target = None       
target = None              

# Action bounds
latest_bounds = None       
# Thrust cmd (ZOH) + bounds at each thrust sample
thr_times, fxs, fys = [], [], []
thr_bounds = []            # per-thrust snapshot (xmin, xmax, ymin, ymax) or None
t0 = None # time zero (first seen message)

def _extract_xy_bounds(ab_msg):
    data = getattr(ab_msg, "data", None)
    if data is None:
        if hasattr(ab_msg, "__iter__"):
            data = list(ab_msg)
        else:
            return None
    data = list(data)
    if len(data) == 6:   # 3x2 -> use first 2 rows
        arr = np.array(data, dtype=float).reshape(3, 2)
    elif len(data) == 4: # 2x2
        arr = np.array(data, dtype=float).reshape(2, 2)
    else:
        return None
    xmin, xmax = float(arr[0, 0]), float(arr[0, 1])
    ymin, ymax = float(arr[1, 0]), float(arr[1, 1])
    return (xmin, xmax, ymin, ymax)

def _is_default_bounds(b, pm=DEFAULT_PM, tol=TOL):
    xmin, xmax, ymin, ymax = b
    x_ok = abs(xmin + pm) <= tol and abs(xmax - pm) <= tol
    y_ok = abs(ymin + pm) <= tol and abs(ymax - pm) <= tol
    return x_ok and y_ok

def _extract_fx_fy_tau(msg):
    if hasattr(msg, "x") and hasattr(msg, "y") and hasattr(msg, "z"):
        return float(msg.x), float(msg.y), float(msg.z)
    for attr in ("vector", "force"):
        sub = getattr(msg, attr, None)
        if sub is not None and all(hasattr(sub, k) for k in ("x", "y", "z")):
            return float(sub.x), float(sub.y), float(sub.z)
    try:
        return float(msg["x"]), float(msg["y"]), float(msg["z"])
    except Exception:
        return None

# Read loop 
while reader.has_next():
    topic, data, t_nsec = reader.read_next()

    if topic == odom_topic:
        msg = deserialize_message(data, Odom)
        if t0 is None:
            t0 = t_nsec
        t = (t_nsec - t0) / 1e9

        times.append(t)
        xs.append(msg.pose.pose.position.x)
        ys.append(msg.pose.pose.position.y)
        vxs.append(msg.twist.twist.linear.x)
        vys.append(msg.twist.twist.linear.y)
        ell_bounds.append(latest_bounds)  # snapshot bounds for this odom

        # ZOH sample of target at odom timestamps
        if latest_target is not None:
            txs.append(latest_target[0])
            tys.append(latest_target[1])
        else:
            txs.append(np.nan)
            tys.append(np.nan)

    elif Target and topic == target_topic:
        tmsg = deserialize_message(data, Target)
        latest_target = (tmsg.pose.pose.position.x, tmsg.pose.pose.position.y)
        target = latest_target  # for XY marker (latest)

    elif ActionBounds and topic == action_bounds_topic:
        ab_msg = deserialize_message(data, ActionBounds)
        b = _extract_xy_bounds(ab_msg)
        if b is not None:
            latest_bounds = b

    elif Thrust and topic == thrust_cmd_topic:
        tmsg = deserialize_message(data, Thrust)
        if t0 is None:
            t0 = t_nsec
        t = (t_nsec - t0) / 1e9

        fx_fy_tau = _extract_fx_fy_tau(tmsg)
        if fx_fy_tau is not None:
            fx, fy, _tau = fx_fy_tau
            thr_times.append(t)
            fxs.append(fx)
            fys.append(fy)
            thr_bounds.append(latest_bounds)  # snapshot bounds for this thrust sample

if not xs:
    raise SystemExit("No odometry data")

# Plot 1: XY trajectory
plt.figure(figsize=(8, 4))
sc = plt.scatter(xs, ys, c=times, cmap="viridis")
# Continuous target trajectory (sampled ZOH at odom timestamps)
if txs and tys:
    tmx = np.array(txs, dtype=float)
    tmy = np.array(tys, dtype=float)
    mask = np.isfinite(tmx) & np.isfinite(tmy)
    if np.any(mask):
        plt.plot(tmx[mask], tmy[mask], label="target traj", linewidth=1)
        # latest target marker
        li = np.where(mask)[0][-1]
        plt.scatter(tmx[li], tmy[li], label="target (latest)", marker="x")

# square [-2,2]
sqx = [-2, 2, 2, -2, -2]
sqy = [-2, -2, 2, 2, -2]
plt.plot(sqx, sqy, linestyle="-.")

# Draw bounds ellipses centered on odom points; axes from bounds, scaled by SCALE
ax = plt.gca()
for x, y, b in zip(xs, ys, ell_bounds):
    if b is None or _is_default_bounds(b):
        continue
    xmin, xmax, ymin, ymax = b
    xmin = np.minimum(0, xmin)
    xmax = np.maximum(0, xmax)
    ymin = np.minimum(0, ymin)
    ymax = np.maximum(0, ymax)
    w = SCALE * abs(xmax - xmin)
    h = SCALE * abs(ymax - ymin)
    e = Ellipse((x, y), width=w, height=h)
    ax.add_patch(e)

plt.axis("equal")
plt.xlabel("$x$ (m)")
plt.ylabel("$y$ (m)")
plt.colorbar(sc, label="Time (s)")
plt.legend()

fig, axes = plt.subplots(2, 1, figsize=(8, 4), sharex=True)

# Plot 2: Position vs time
axes[0].plot(times, xs, label="$x$")
axes[0].plot(times, ys, label="$y$")
axes[0].axhline(2, linestyle=":", label=r"$\mathbf{b}$")
axes[0].axhline(-2, linestyle=":")
axes[0].set_ylabel("Position (m)")
axes[0].legend()

# Plot 3: Velocity vs time
axes[1].plot(times, vxs, label="$v_x$")
axes[1].plot(times, vys, label="$v_y$")
axes[1].set_xlabel("Time (s)")
axes[1].set_ylabel("Velocity (m/s)")
axes[1].legend()

# Plot 4: Thrust cmds (2-row) with bounds (ZOH)
if thr_times:
    xmins, xmaxs, ymins, ymaxs = [], [], [], []
    for b in thr_bounds:
        if b is None:
            xmins.append(np.nan); xmaxs.append(np.nan)
            ymins.append(np.nan); ymaxs.append(np.nan)
        else:
            xmin, xmax, ymin, ymax = b
            xmins.append(xmin); xmaxs.append(xmax)
            ymins.append(ymin); ymaxs.append(ymax)

    fig, (ax_fx, ax_fy) = plt.subplots(2, 1, sharex=True, figsize=(8, 4))
    ax_fx.step(thr_times, fxs, where="post", label="$f_x$", linewidth=1)
    ax_fx.step(thr_times, xmins, where="post", label=r"$f_{x,\mathrm{min}}$", linewidth=0.9, linestyle=":")
    ax_fx.step(thr_times, xmaxs, where="post", label=r"$f_{x,\mathrm{max}}$", linewidth=0.9, linestyle=":")
    ax_fx.set_ylabel("$f_x$ (N)")
    ax_fx.legend()

    ax_fy.step(thr_times, fys, where="post", label="$f_y$", linewidth=1)
    ax_fy.step(thr_times, ymins, where="post", label=r"$f_{y,\mathrm{min}}$", linewidth=0.9, linestyle=":")
    ax_fy.step(thr_times, ymaxs, where="post", label=r"$f_{y,\mathrm{max}}$", linewidth=0.9, linestyle=":")
    ax_fy.set_xlabel("Time (s)")
    ax_fy.set_ylabel("$f_y$ (N)")
    ax_fy.legend()

# Plot 5: Target vs Odom and Error (3-row: x, y, L2)
if txs and np.isfinite(np.array(txs)).any():
    xs_arr = np.array(xs); ys_arr = np.array(ys)
    txs_arr = np.array(txs); tys_arr = np.array(tys)

    ex = txs_arr - xs_arr
    ey = tys_arr - ys_arr
    eL2 = np.sqrt(ex**2 + ey**2)

    fig, (ax_x, ax_y, ax_l2) = plt.subplots(3, 1, sharex=True, figsize=(9, 4))
    ax_x.plot(times, ex)
    ax_x.axhline(0, linestyle=":")
    ax_x.set_ylabel(r"$\delta x$ (m)")
    ax_x.legend()

    ax_y.plot(times, ey)
    ax_y.axhline(0, linestyle=":")
    ax_y.set_ylabel(r"$\delta y$ (m)")
    ax_y.legend()

    ax_l2.plot(times, eL2)
    ax_l2.axhline(0, linestyle=":")
    ax_l2.set_xlabel("Time (s)")
    ax_l2.set_ylabel(r"$\||\delta r\||$ (m)")
    ax_l2.legend()

# ---------------- Save figs ----------------
out_dir = os.path.dirname(bag) or "."
for i, num in enumerate(sorted(plt.get_fignums()), start=1):
    fig = plt.figure(num)
    fig.savefig(os.path.join(out_dir, f"fig_{i}.png"))
plt.close("all")
