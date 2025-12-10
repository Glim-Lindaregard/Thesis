#include "slider_controller.hpp"

#include <ignition/math/Vector3.hh>

namespace gazebo_plugins
{

SliderController::SliderController() : gazebo::ModelPlugin() {
}

SliderController::~SliderController() {
}

void SliderController::Load(gazebo::physics::ModelPtr model, sdf::ElementPtr sdf) {
    // Get ROS node
    this->ros_node = gazebo_ros::Node::Get(sdf);

    // Get QoS profile
    const gazebo_ros::QoS & qos = this->ros_node->get_qos();

    // Check for manual force
    if(sdf->HasElement("force")) {
        this->force = sdf->GetElement("force")->Get<double>();
    } else {
        RCLCPP_WARN(rclcpp::get_logger("slider_controller"), "No thruster force specified, using default value of %f", this->DEFAULT_FORCE);
        this->force = this->DEFAULT_FORCE;
    }

    // Check for manual torque
    if(sdf->HasElement("torque")) {
        this->torque = sdf->GetElement("torque")->Get<double>();
    } else {
        RCLCPP_WARN(rclcpp::get_logger("slider_controller"), "No thruster torque specified, using default value of %f", this->DEFAULT_TORQUE);
        this->torque = this->DEFAULT_TORQUE;
    }

    // Get the base link
    this->cog_link = model->GetLinks()[0];
    if(!this->cog_link){
        RCLCPP_ERROR(rclcpp::get_logger("slider_controller"), "Link not found! Exiting...");
        return;
    } else {
        RCLCPP_DEBUG(rclcpp::get_logger("slider_controller"), "Link <%s> found!", this->cog_link->GetName().c_str());
    }

    // Create subscriber
    this->sub = this->ros_node->create_subscription<std_msgs::msg::UInt8MultiArray>(
        "slider_controller", qos.get_subscription_qos("controller", rclcpp::SystemDefaultsQoS()),
        std::bind(&SliderController::OnMsg, this, std::placeholders::_1));

    // Callback on every world update
    this->update_connection = gazebo::event::Events::ConnectWorldUpdateBegin(
        std::bind(&SliderController::OnUpdate, this));

    RCLCPP_INFO(rclcpp::get_logger("slider_controller"), "Slider Controller Plugin loaded on %s!\n", model->GetName().c_str());
}

void SliderController::OnMsg(const std_msgs::msg::UInt8MultiArray::SharedPtr msg) {
    // Copy forces from ROS message
    for (int i = 0; i < 8; i++) {
        this->forces[i] = msg->data[i];
    }

    // Print forces
    RCLCPP_DEBUG(rclcpp::get_logger("slider_controller"), "Forces: %d %d %d %d %d %d %d %d", this->forces[0], this->forces[1], this->forces[2], this->forces[3], this->forces[4], this->forces[5], this->forces[6], this->forces[7]);
}

void SliderController::OnUpdate() {
    // Apply forces
    for (int i = 0; i < 8; i++) {
        if(this->forces[i]){
            this->cog_link->AddRelativeForce(
                this->force * this->force_vectors[i]
            );
            
            this->cog_link->AddRelativeTorque(
                this->torque * this->torque_vectors[i]
            );
        }
    }
}

GZ_REGISTER_MODEL_PLUGIN(SliderController)
} // namespace gazebo