from isaaclab.utils import configclass
from isaaclab.scene import InteractiveSceneCfg
from isaaclab.sensors import CameraCfg, ContactSensorCfg
from isaaclab.terrains import TerrainImporterCfg
from isaaclab.assets import AssetBaseCfg, RigidObjectCfg
import isaaclab.sim as sim_utils
from isaaclab.sensors import ImuCfg
from pathlib import Path

import os

# Local paths to usda models
current_dir_path = os.path.dirname(os.path.abspath(__file__))
slider_usd_path = "../../assets/body.usda"
slider_usd_path = os.path.join(current_dir_path, slider_usd_path)


@configclass
class SliderSceneCfg(InteractiveSceneCfg):
    """
    Defines a 4x4 box arena centered at (0,0) with the slider robot 
    and target point represented by as sphere 
    
    """
    
    # Ground plane and light source 
    terrain = TerrainImporterCfg(prim_path="/World/ground", 
                                    terrain_type="plane", 
                                    debug_vis=False, 
                                    visual_material=sim_utils.PreviewSurfaceCfg(diffuse_color=(0.0, 0.0, 0.0)),
                                    collision_group=-1
                                    )
    
    light = AssetBaseCfg(
        prim_path="/World/light",
        spawn=sim_utils.DomeLightCfg(
            intensity=2000.0, 
            color=(0.8, 0.8, 0.8)
        )
    )

    target_sphere: RigidObjectCfg = RigidObjectCfg(
        prim_path="{ENV_REGEX_NS}/target_sphere",
        spawn=sim_utils.SphereCfg(
            radius=0.1, 
            rigid_props=sim_utils.RigidBodyPropertiesCfg(kinematic_enabled=False),
            collision_props=sim_utils.CollisionPropertiesCfg(
                collision_enabled=False, 
                rest_offset=0.0,
                contact_offset=0.1
                ), 
            visual_material=sim_utils.PreviewSurfaceCfg(diffuse_color=(0.5, 0.0, 0.0), 
                                                        metallic=0.0)
        ), 
        init_state=RigidObjectCfg.InitialStateCfg(
            pos=(0.0, 0.0, 0.1),
            rot=(1.0, 0.0, 0.0, 0.0)
        ),
        collision_group=0,
    )

    slider: RigidObjectCfg = RigidObjectCfg(
        prim_path="{ENV_REGEX_NS}/slider", 
        spawn=sim_utils.UsdFileCfg(
            usd_path=slider_usd_path, 
            visual_material=sim_utils.PreviewSurfaceCfg(
                diffuse_color=(0.1, 0.1, 0.1), 
                metallic=0.2
            ), 
            mass_props=sim_utils.MassPropertiesCfg(mass=4.528),
            activate_contact_sensors=True, 
            collision_props=sim_utils.CollisionPropertiesCfg(
                collision_enabled=True, 
                rest_offset=0.0,
                contact_offset=0.1
            )
        ), 
        init_state=RigidObjectCfg.InitialStateCfg(
            pos=(0.0, 0.0, 0.1), 
            rot=(1.0, 0.0, 0.0, 0.0)
        ),
        collision_group=0,
        debug_vis=True
    )
    
    # Make arena walls
    
    wall_n = RigidObjectCfg(
        prim_path="{ENV_REGEX_NS}/wall_n",
        spawn=sim_utils.CuboidCfg(
            size=(5.2, 0.2, 1),
            rigid_props=sim_utils.RigidBodyPropertiesCfg(kinematic_enabled=True),
            mass_props=sim_utils.MassPropertiesCfg(mass=1e9),
            collision_props=sim_utils.CollisionPropertiesCfg(
                collision_enabled=True,
                rest_offset=0.0, 
                contact_offset=0.01
            ),
            visual_material=sim_utils.PreviewSurfaceCfg(
                diffuse_color=(0.0, 0.0, 1.0), 
                metallic=0.2
            ),
            physics_material=sim_utils.RigidBodyMaterialCfg(),
        ),
        init_state=RigidObjectCfg.InitialStateCfg(pos=(0.0, 2.5, 0.0), rot=(1.0, 0.0, 0.0, 0.0)),
    )

    wall_s = RigidObjectCfg(
        prim_path="{ENV_REGEX_NS}/wall_s",
        spawn=sim_utils.CuboidCfg(
            size=(5.2, 0.2, 1),
            rigid_props=sim_utils.RigidBodyPropertiesCfg(kinematic_enabled=True),
            mass_props=sim_utils.MassPropertiesCfg(mass=1e9),
            collision_props=sim_utils.CollisionPropertiesCfg(
                collision_enabled=True, 
                rest_offset=0.0, 
                contact_offset=0.01
            ),
            visual_material=sim_utils.PreviewSurfaceCfg(
                diffuse_color=(0.0, 0.0, 1.0), 
                metallic=0.2
            ),
            physics_material=sim_utils.RigidBodyMaterialCfg(),
        ),
        init_state=RigidObjectCfg.InitialStateCfg(
            pos=(0.0, -2.5, 0.0),
            rot=(1.0, 0.0, 0.0, 0.0)
        ),
    )

    wall_e = RigidObjectCfg(
        prim_path="{ENV_REGEX_NS}/wall_e",
        spawn=sim_utils.CuboidCfg(
            size=(0.2, 5.2, 1),
            rigid_props=sim_utils.RigidBodyPropertiesCfg(kinematic_enabled=True),
            mass_props=sim_utils.MassPropertiesCfg(mass=1e9),
            collision_props=sim_utils.CollisionPropertiesCfg(
                collision_enabled=True,
                rest_offset=0.0,
                contact_offset=0.01
            ),
            visual_material=sim_utils.PreviewSurfaceCfg(
                diffuse_color=(0.0, 0.0, 1.0),
                metallic=0.2
            ),
            physics_material=sim_utils.RigidBodyMaterialCfg(),
        ),
        init_state=RigidObjectCfg.InitialStateCfg(
            pos=(2.5, 0.0, 0.0), 
            rot=(1.0, 0.0, 0.0, 0.0)
        ),
    )

    wall_w = RigidObjectCfg(
        prim_path="{ENV_REGEX_NS}/wall_w",
        spawn=sim_utils.CuboidCfg(
            size=(0.2, 5.2, 1),
            rigid_props=sim_utils.RigidBodyPropertiesCfg(kinematic_enabled=True),
            mass_props=sim_utils.MassPropertiesCfg(mass=1e9),
            collision_props=sim_utils.CollisionPropertiesCfg(
                collision_enabled=True, 
                rest_offset=0.0,
                contact_offset=0.01
            ),
            visual_material=sim_utils.PreviewSurfaceCfg(
                diffuse_color=(0.0, 0.0, 1.0),
                metallic=0.2
            ),
            physics_material=sim_utils.RigidBodyMaterialCfg(),
        ),
        init_state=RigidObjectCfg.InitialStateCfg(
            pos=(-2.5, 0.0, 0.0),
            rot=(1.0, 0.0, 0.0, 0.0)
        ),
    )



    # Make sensors
    # camera = CameraCfg(
    #     prim_path="{ENV_REGEX_NS}/slider/body/front_cam",
    #     update_period=0.1,
    #     height=445,
    #     width=445,
    #     data_types=["rgb", "distance_to_image_plane"],
    #     spawn=sim_utils.PinholeCameraCfg(
    #         focal_length=24.0, 
    #         focus_distance=400.0,
    #         horizontal_aperture=20.955, 
    #         clipping_range=(0.1, 1.0e5)
    #     ),
    #     offset=CameraCfg.OffsetCfg(pos=(0.12, 0.0, 0.019),
    #                                rot=(0.5, -0.5, 0.5, -0.5),
    #                                convention="ros"),
    #     debug_vis=False
    # )

    contact_forces = ContactSensorCfg(
        prim_path="{ENV_REGEX_NS}/slider/body", 
        update_period=0.0,
        history_length=6,
        debug_vis=False
    )
 
    imu = ImuCfg(
        prim_path="{ENV_REGEX_NS}/slider/body", 
        gravity_bias=(0.0, 0.0, 9.81),
        debug_vis=True
    )
