#ifndef GAZEBO_SLIDER_CONTROLLER_HPP_
#define GAZEBO_SLIDER_CONTROLLER_HPP_

#include <gazebo/common/Events.hh>
#include <gazebo/physics/Link.hh>
#include <gazebo/physics/Model.hh>
#include <gazebo/physics/physics.hh>
#include <rclcpp/rclcpp.hpp>
#include <gazebo_ros/node.hpp>
#include <gazebo_ros/utils.hpp>
#include <gazebo_ros/conversions/geometry_msgs.hpp>
#include <gazebo_ros/conversions/builtin_interfaces.hpp>
#include "std_msgs/msg/u_int8_multi_array.hpp"
#include <gazebo/common/Plugin.hh>

using namespace std;

namespace gazebo_plugins
{

class SliderController : public gazebo::ModelPlugin {
    public:
        // Constructor
        SliderController();

        // Destructor
        virtual ~SliderController();

    protected:
        // On loading the plugin
        void Load(gazebo::physics::ModelPtr model, sdf::ElementPtr sdf) override;
        
        // On world update
        virtual void OnUpdate();

    private:
        // On message from ROS
        void OnMsg(const std_msgs::msg::UInt8MultiArray::SharedPtr msg);

        // Pointer to the link where the COG is
        gazebo::physics::LinkPtr cog_link;

        // gazebo::physics::LinkPtr links[8];

        // Pointer to the GazeboROS node
        gazebo_ros::Node::SharedPtr ros_node;

        // Subscriber
        rclcpp::Subscription<std_msgs::msg::UInt8MultiArray>::SharedPtr sub;

        // Container for the forces to be applied from T11 to T42
        uint8_t forces[8];

        // Pointer to the update event connection
        gazebo::event::ConnectionPtr update_connection;

        // Thruster Force & Torque
        double force;
        double torque;
        const double DEFAULT_FORCE = 0.7;
        const double DEFAULT_TORQUE = 0.0098;

        // Force direction vectors
        const ignition::math::Vector3<double> force_vectors[8] = {
            ignition::math::Vector3<double>(-1.0f,  0.0f, 0.0f),
            ignition::math::Vector3<double>( 0.0f,  1.0f, 0.0f),
            ignition::math::Vector3<double>( 1.0f,  0.0f, 0.0f),
            ignition::math::Vector3<double>( 0.0f,  1.0f, 0.0f),
            ignition::math::Vector3<double>( 1.0f,  0.0f, 0.0f),
            ignition::math::Vector3<double>( 0.0f, -1.0f, 0.0f),
            ignition::math::Vector3<double>(-1.0f,  0.0f, 0.0f),
            ignition::math::Vector3<double>( 0.0f, -1.0f, 0.0f),
        };

        // Torque direction vectors
        const ignition::math::Vector3<double> torque_vectors[8] = {
            ignition::math::Vector3<double>(0.0f, 0.0f, -1.0f),
            ignition::math::Vector3<double>(0.0f, 0.0f,  1.0f),
            ignition::math::Vector3<double>(0.0f, 0.0f,  1.0f),
            ignition::math::Vector3<double>(0.0f, 0.0f, -1.0f),
            ignition::math::Vector3<double>(0.0f, 0.0f, -1.0f),
            ignition::math::Vector3<double>(0.0f, 0.0f,  1.0f),
            ignition::math::Vector3<double>(0.0f, 0.0f,  1.0f),
            ignition::math::Vector3<double>(0.0f, 0.0f, -1.0f),
        };
};

}

#endif