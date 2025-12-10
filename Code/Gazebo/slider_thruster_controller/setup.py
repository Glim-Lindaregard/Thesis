from setuptools import setup
import os 
from glob import glob 


package_name = 'slider_thruster_controller'
tsl_optimizer = 'tsl_optimizer'


setup(
    name=package_name,
    version='0.0.0',
    packages=[package_name, tsl_optimizer],
    # install_requires=[
    #     'setuptools',
    #     'quadprog', 
    #     'cvxopt', 
    #     'numpy'
    # ],
    data_files=[
        ('share/ament_index/resource_index/packages',
            ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
        (os.path.join('share', package_name, 'launch'), glob(os.path.join('launch', '*.launch')))      
    ],
    install_requires=['setuptools'],
    zip_safe=True,
    maintainer='bibeto',
    maintainer_email='bibeto@todo.todo',
    description='TODO: Package description',
    license='TODO: License declaration',
    tests_require=['pytest'],
    entry_points={
        'console_scripts': [
            'pwm_publisher = slider_thruster_controller.pwm_publisher:main',
            'pid_controller = slider_thruster_controller.pid_controller:main',
            'pid_parameters = slider_thruster_controller.pid_parameters:main',
            'solver = ' + tsl_optimizer + '.solver:main', 
            'parameters = ' + tsl_optimizer + '.parameters:main'
        ],
    },
)
