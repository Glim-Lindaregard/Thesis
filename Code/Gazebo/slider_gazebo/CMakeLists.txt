cmake_minimum_required(VERSION 3.5)
project(slider_gazebo)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclpy REQUIRED)
find_package(gazebo_ros REQUIRED)


# plugin dependencies 
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(gazebo_dev REQUIRED)
find_package(gazebo_msgs REQUIRED)
find_package(gazebo_ros REQUIRED)


link_directories(${gazebo_dev_LIBRARY_DIRS})

include_directories(include
  ${rclcpp_INCLUDE_DIRS}
  ${std_msgs_INCLUDE_DIRS}
  ${gazebo_ros_INCLUDE_DIRS}
)

find_package(gazebo REQUIRED)
include_directories(${GAZEBO_INCLUDE_DIRS})
link_directories(${GAZEBO_LIBRARY_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GAZEBO_CXX_FLAGS}")

add_library(slider_controller SHARED plugin/slider_controller.cpp)
target_include_directories(slider_controller PUBLIC include plugin)
ament_target_dependencies(slider_controller 
  "gazebo_dev" 
  "gazebo_ros"
  "rclcpp" 
  "geometry_msgs" 
  "gazebo_msgs" 
  "std_msgs"
) 
ament_export_libraries(slider_controller)


install(TARGETS
    slider_controller
  DESTINATION lib
)
# plugin specific CMake commands end here


if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  ament_lint_auto_find_test_dependencies()
endif()

install(
  DIRECTORY 
    launch 
    worlds
  DESTINATION
    share/${PROJECT_NAME}/
)

install(
  PROGRAMS
  launch/urdf_publisher.py
  DESTINATION lib/${PROJECT_NAME}
)


ament_package()