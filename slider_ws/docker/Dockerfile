FROM osrf/ros:humble-desktop

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get upgrade -y

RUN apt-get update && apt-get install -y ros-humble-joint-state-publisher
RUN apt-get update && apt-get install -y ros-humble-tf-transformations

RUN apt-get update && apt-get install -y \
    gazebo \
    ros-humble-gazebo-ros-pkgs \
    ros-humble-xacro

RUN apt-get update && apt-get install -y \
    iproute2 \
    iputils-ping \
    net-tools \
    dnsutils \
    curl \
    wget 


RUN apt-get install -y python3-pip

RUN pip install setuptools quadprog cvxopt numpy stable-baselines3 casadi \
    && python3 -m pip install --upgrade "numpy<1.24"


RUN apt-get update && apt-get install -y \
    ros-humble-plotjuggler-ros \
    btop

RUN apt-get install -y wget git p7zip-full




RUN echo 'source /opt/ros/humble/setup.bash' >> /etc/bash.bashrc && \
    echo '[ -f /root/slider_ws/install/setup.bash ] && source /root/slider_ws/install/setup.bash' >> /etc/bash.bashrc

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
CMD ["bash"]

