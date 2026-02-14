import os
from pathlib import Path

import bpy
from bpy.props import FloatProperty, PointerProperty, StringProperty
from bpy.types import Operator, Panel, PropertyGroup
from mathutils import Vector


def scene_get(key, default):
    return lambda settings: settings.id_data.get(key, default)


def scene_set(key):
    return lambda settings, value: settings.id_data.__setitem__(key, value)


def camera_get(key, default):
    def get(settings):
        camera = settings.id_data.camera
        return camera.data.get(key, default) if camera else default
    return get


def camera_set(key):
    def set_value(settings, value):
        camera = settings.id_data.camera
        if camera:
            camera.data[key] = value
    return set_value


class RTSettings(PropertyGroup):
    export_path: StringProperty(name="GLB", subtype="FILE_PATH")
    hdri: StringProperty(
        name="HDRI", subtype="FILE_PATH",
        get=scene_get("raytracer_hdri", ""),
        set=scene_set("raytracer_hdri"))
    hdri_rotation: FloatProperty(
        name="Rotation", subtype="ANGLE",
        get=scene_get("raytracer_hdri_rotation", 0.0),
        set=scene_set("raytracer_hdri_rotation"))
    hdri_strength: FloatProperty(
        name="Strength", min=0.0,
        get=scene_get("raytracer_hdri_strength", 1.0),
        set=scene_set("raytracer_hdri_strength"))
    exposure: FloatProperty(
        name="Exposure",
        get=scene_get("raytracer_exposure", 0.0),
        set=scene_set("raytracer_exposure"))
    aperture: FloatProperty(
        name="Lens radius", min=0.0, unit="LENGTH",
        get=camera_get("raytracer_aperture", 0.0),
        set=camera_set("raytracer_aperture"))
    focus_distance: FloatProperty(
        name="Focus", min=0.001, unit="LENGTH",
        get=camera_get("raytracer_focus_distance", 10.0),
        set=camera_set("raytracer_focus_distance"))


def target(context):
    obj = context.active_object
    if not obj:
        return None
    if obj.type == "MESH":
        return sum((obj.matrix_world @ Vector(corner) for corner in obj.bound_box),
                   Vector()) / 8.0
    return obj.matrix_world.translation


class RT_OT_focus_selected(Operator):
    bl_idname = "raytracer.focus_selected"
    bl_label = "Focus selected"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return context.scene.camera and target(context) is not None

    def execute(self, context):
        camera = context.scene.camera
        distance = (target(context) - camera.matrix_world.translation).length
        context.scene.raytracer_tools.focus_distance = distance
        camera.data.dof.focus_distance = distance
        return {"FINISHED"}


class RT_OT_aim_selected(Operator):
    bl_idname = "raytracer.aim_selected"
    bl_label = "Aim and focus"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return context.scene.camera and target(context) is not None

    def execute(self, context):
        camera = context.scene.camera
        point = target(context)
        camera.rotation_euler = (
            point - camera.matrix_world.translation).to_track_quat("-Z", "Y").to_euler()
        bpy.ops.raytracer.focus_selected()
        return {"FINISHED"}


class RT_OT_export(Operator):
    bl_idname = "raytracer.export"
    bl_label = "Export GLB"

    def execute(self, context):
        settings = context.scene.raytracer_tools
        source = settings.export_path or bpy.data.filepath
        if not source:
            self.report({"ERROR"}, "Save the blend file or choose an export path")
            return {"CANCELLED"}
        path = Path(bpy.path.abspath(source))
        path = path.with_suffix(".glb")
        path.parent.mkdir(parents=True, exist_ok=True)

        hdri = context.scene.get("raytracer_hdri", "")
        if hdri:
            if hdri.startswith("//"):
                source = Path(bpy.path.abspath(hdri))
            elif Path(hdri).is_absolute():
                source = Path(hdri)
            else:
                source = Path(bpy.data.filepath).parent / hdri
            context.scene["raytracer_hdri"] = os.path.relpath(
                source, path.parent).replace(os.sep, "/")

        bpy.ops.export_scene.gltf(
            filepath=str(path), export_format="GLB", export_cameras=True,
            export_lights=True, export_extras=True, export_apply=True)
        self.report({"INFO"}, f"Exported {path.name}")
        return {"FINISHED"}


class RT_PT_tools(Panel):
    bl_label = "Raytracer"
    bl_idname = "RT_PT_tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Raytracer"

    def draw(self, context):
        layout = self.layout
        settings = context.scene.raytracer_tools
        layout.prop(settings, "export_path")
        layout.prop(settings, "hdri")
        row = layout.row(align=True)
        row.prop(settings, "hdri_rotation")
        row.prop(settings, "hdri_strength")
        layout.prop(settings, "exposure")
        layout.separator()
        column = layout.column()
        column.enabled = context.scene.camera is not None
        column.prop(settings, "aperture")
        column.prop(settings, "focus_distance")
        row = column.row(align=True)
        row.operator("raytracer.focus_selected")
        row.operator("raytracer.aim_selected")
        layout.separator()
        layout.operator("raytracer.export", icon="EXPORT")


classes = (
    RTSettings,
    RT_OT_focus_selected,
    RT_OT_aim_selected,
    RT_OT_export,
    RT_PT_tools,
)


def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.raytracer_tools = PointerProperty(type=RTSettings)


def unregister():
    del bpy.types.Scene.raytracer_tools
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
