"""Render the single-fibre Re=200 HTG surface/PD points view at 49 s.

Run with:

    pvpython generate_2.py

By default this writes one PNG into ../../raster.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from paraview.simple import (
    ColorBy,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetLayout,
    OpenDataFile,
    Render,
    SaveScreenshot,
    SetActiveSource,
    Show,
    _DisableFirstRenderCameraReset,
)


ROOT = Path(__file__).resolve().parent
CASE_DIR = ROOT / "case"
DEFAULT_OUTPUT_DIR = ROOT.parent.parent / "raster"
DEFAULT_OUTPUT_NAME = f"{ROOT.name}_t49_surface_edges.png"

HTG_SERIES = CASE_DIR / "htg.vtkhdf.series"
PD_SERIES = CASE_DIR / "pd.vtkhdf.series"

TARGET_TIME = 49.0
IMAGE_RESOLUTION = (1600, 900)
VIEW_SIZE = (1600, 900)
BACKGROUND = (1.0, 1.0, 1.0)

HTG_U_FIELD = "u"
HTG_U_COMPONENT = "Magnitude"
HTG_EDGE_COLOR = [0.0, 0.0, 0.0]
HTG_EDGE_LINE_WIDTH = 0.05

PD_POINT_COLOR = [1.0, 0.0, 0.0]
PD_POINT_SIZE = 10.0


def parse_args():
    parser = argparse.ArgumentParser(
        description="Render the 49 s HTG surface-with-edges view.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Directory for the generated PNG file "
            f"(default: {DEFAULT_OUTPUT_DIR})"
        ),
    )
    parser.add_argument(
        "--output-name",
        default=DEFAULT_OUTPUT_NAME,
        help=f"Generated PNG filename (default: {DEFAULT_OUTPUT_NAME}).",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=TARGET_TIME,
        help=f"Requested render time; nearest available time is used (default: {TARGET_TIME}).",
    )
    parser.add_argument(
        "--resolution",
        type=int,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=IMAGE_RESOLUTION,
        help=(
            "Saved PNG resolution in pixels "
            f"(default: {IMAGE_RESOLUTION[0]} {IMAGE_RESOLUTION[1]})."
        ),
    )
    parser.add_argument(
        "--view-size",
        type=int,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=VIEW_SIZE,
        help=(
            "Live ParaView view size used before tiled screenshot export "
            f"(default: {VIEW_SIZE[0]} {VIEW_SIZE[1]})."
        ),
    )
    return parser.parse_args()


def create_view(view_size):
    view = GetActiveViewOrCreate("RenderView")

    layout = GetLayout()
    layout.SetSize(*view_size)

    view.ViewSize = view_size
    view.Background = BACKGROUND
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.Set(
        InteractionMode="2D",
        CameraPosition=[1.03767530860879, 0.08054220125218091, 43.71281292110204],
        CameraFocalPoint=[1.03767530860879, 0.08054220125218091, 0.0],
        CameraViewUp=[-1.0, 0.0, 0.0],
        CameraParallelScale=1.6817106776966735,
    )
    return view


def build_pipeline(view):
    htg = OpenDataFile(str(HTG_SERIES))
    pd = OpenDataFile(str(PD_SERIES))

    htg_display = Show(htg, view, "HyperTreeGridRepresentation")
    htg_display.SetRepresentationType("Surface With Edges")
    ColorBy(htg_display, ("CELLS", HTG_U_FIELD, HTG_U_COMPONENT))
    htg_display.EdgeColor = HTG_EDGE_COLOR
    htg_display.LineWidth = HTG_EDGE_LINE_WIDTH
    htg_display.SetScalarBarVisibility(view, False)

    pd_display = Show(pd, view, "GeometryRepresentation")
    pd_display.SetRepresentationType("Points")
    ColorBy(pd_display, ("POINTS", None))
    pd_display.AmbientColor = PD_POINT_COLOR
    pd_display.DiffuseColor = PD_POINT_COLOR
    pd_display.PointSize = PD_POINT_SIZE
    pd_display.RenderPointsAsSpheres = 0

    SetActiveSource(htg)
    return htg, pd, htg_display


def available_times(source):
    source.UpdatePipeline()
    return [float(time) for time in source.TimestepValues]


def nearest_time(times, target):
    return min(times, key=lambda time: abs(time - target))


def render_time(view, sources, htg_display, time, output_path, image_resolution):
    for source in sources:
        source.UpdatePipeline(time)

    htg_display.RescaleTransferFunctionToDataRange(False, False)

    scene = GetAnimationScene()
    scene.AnimationTime = time
    view.ViewTime = time
    Render(view)

    SaveScreenshot(
        str(output_path),
        view,
        ImageResolution=image_resolution,
        TransparentBackground=1,
    )


def main():
    args = parse_args()

    _DisableFirstRenderCameraReset()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    view = create_view(args.view_size)
    htg, pd, htg_display = build_pipeline(view)
    time = nearest_time(available_times(htg), args.time)

    output_path = args.output_dir / args.output_name
    print(
        f"Rendering t = {time:.12g} at {args.resolution[0]}x{args.resolution[1]} "
        f"-> {output_path}"
    )
    render_time(
        view,
        (htg, pd),
        htg_display,
        time,
        output_path,
        args.resolution,
    )


if __name__ == "__main__":
    main()
