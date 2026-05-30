"""Render selected time steps from the single-fibre Re=200 ParaView series.

Run with:

    pvpython generate.py --output-dir ../../raster

The script currently renders only the final stored time step.  Add more
entries to TIME_GROUPS to write additional sets of frames.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from paraview.simple import (
    ColorBy,
    Contour,
    GetAnimationScene,
    GetActiveViewOrCreate,
    GetColorTransferFunction,
    GetLayout,
    GetOpacityTransferFunction,
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
DEFAULT_OUTPUT_PREFIX = ROOT.name

HTG_SERIES = CASE_DIR / "htg.vtkhdf.series"
PD_SERIES = CASE_DIR / "pd.vtkhdf.series"

IMAGE_RESOLUTION = (425, 4000)
BACKGROUND = (1.0, 1.0, 1.0)
CONTOUR_COLOR = [0.0, 0.0, 0.0]
PD_POINT_COLOR = [1.0, 0.0, 0.0]
PD_POINT_SIZE = 6.0

OMEGA_FIELD = "omega"
OMEGA_COLOR_RANGE = (-4.0, 4.0)
NEGATIVE_OMEGA_LEVELS = [-4.0, -3.25, -2.5, -1.75, -1.0, -0.25]
POSITIVE_OMEGA_LEVELS = [0.25, 1.0, 1.75, 2.5, 3.25, 4.0]

# Available selectors:
#   {"indices": [-1]}                 select by Python indices
#   {"times": [10.0, 20.0, 30.0]}     select nearest available times
#   {"range": (20.0, 30.0, 2.5)}      select nearest times over a range

TIME_GROUPS = {
    "early": {"times": [9.2,10.0,10.8,11.6]},
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Render selected time steps from a ParaView time series.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Directory for generated PNG files "
            f"(default: {DEFAULT_OUTPUT_DIR})"
        ),
    )
    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        help=(
            "Filename prefix for generated PNG files. Use an empty string to "
            f"disable it (default: {DEFAULT_OUTPUT_PREFIX})."
        ),
    )
    return parser.parse_args()


def create_view():
    view = GetActiveViewOrCreate("RenderView")

    layout = GetLayout()
    layout.SetSize(*IMAGE_RESOLUTION)

    view.ViewSize = IMAGE_RESOLUTION
    view.Background = BACKGROUND
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.Set(
        InteractionMode="2D",
        CameraPosition=[4.0, 0.0, 43.71281292110204],
        CameraFocalPoint=[4.0, 0.0, 0.0],
        CameraViewUp=[-1.0, 0.0, 0.0],
        CameraParallelScale=5.277928512317948,
    )
    return view


def build_pipeline(view):
    htg = OpenDataFile(str(HTG_SERIES))
    pd = OpenDataFile(str(PD_SERIES))

    # htg_display = Show(htg, view, "HyperTreeGridRepresentation")
    # htg_display.SetRepresentationType("HTG Surface")
    # ColorBy(htg_display, ("CELLS", OMEGA_FIELD))
    # htg_display.SetScalarBarVisibility(view, False)
    # configure_omega_colormap()

    negative_contour = Contour(registrationName="Negative Omega", Input=htg)
    negative_contour.ContourBy = ["CELLS", OMEGA_FIELD]
    negative_contour.Isosurfaces = NEGATIVE_OMEGA_LEVELS
    negative_display = Show(negative_contour, view, "GeometryRepresentation")
    negative_display.Representation = "Surface"
    negative_display.DiffuseColor = CONTOUR_COLOR
    negative_display.LineWidth = 1.5

    positive_contour = Contour(registrationName="Positive Omega", Input=htg)
    positive_contour.ContourBy = ["CELLS", OMEGA_FIELD]
    positive_contour.Isosurfaces = POSITIVE_OMEGA_LEVELS
    positive_display = Show(positive_contour, view, "GeometryRepresentation")
    positive_display.Representation = "Surface"
    positive_display.DiffuseColor = CONTOUR_COLOR
    positive_display.LineWidth = 1.5

    pd_display = Show(pd, view, "GeometryRepresentation")
    pd_display.Representation = "Points"
    pd_display.AmbientColor = PD_POINT_COLOR
    pd_display.DiffuseColor = PD_POINT_COLOR
    pd_display.PointSize = PD_POINT_SIZE
    pd_display.RenderPointsAsSpheres = 0

    SetActiveSource(htg)
    return htg, pd


def configure_omega_colormap():
    omega_lut = GetColorTransferFunction(OMEGA_FIELD)
    omega_lut.RescaleTransferFunction(*OMEGA_COLOR_RANGE)

    omega_pwf = GetOpacityTransferFunction(OMEGA_FIELD)
    omega_pwf.RescaleTransferFunction(*OMEGA_COLOR_RANGE)


def available_times(source):
    source.UpdatePipeline()
    return [float(time) for time in source.TimestepValues]


def nearest_time(times, target):
    return min(times, key=lambda time: abs(time - target))


def select_time_group(times, selection):
    selected = []

    for index in selection.get("indices", []):
        selected.append(times[index])

    for target in selection.get("times", []):
        selected.append(nearest_time(times, target))

    if "range" in selection:
        start, stop, step = selection["range"]
        target = start
        while target <= stop + 0.5 * step:
            selected.append(nearest_time(times, target))
            target += step

    return unique_in_order(selected)


def unique_in_order(values):
    unique = []
    for value in values:
        if value not in unique:
            unique.append(value)
    return unique


def selected_times(times):
    for group_name, selection in TIME_GROUPS.items():
        for time in select_time_group(times, selection):
            yield group_name, time


def time_label(time):
    label = f"{time:.6f}".rstrip("0").rstrip(".")
    return label.replace("-", "m").replace(".", "p")


def output_filename(prefix, group_name, time):
    stem = f"{group_name}_t{time_label(time)}"
    if prefix:
        stem = f"{prefix}_{stem}"
    return f"{stem}.png"


def render_time(view, sources, time, output_path):
    for source in sources:
        source.UpdatePipeline(time)

    scene = GetAnimationScene()
    scene.AnimationTime = time
    view.ViewTime = time
    Render(view)

    SaveScreenshot(
        str(output_path),
        view,
        ImageResolution=IMAGE_RESOLUTION,
        TransparentBackground=1,
    )


def main():
    args = parse_args()

    _DisableFirstRenderCameraReset()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    view = create_view()
    htg, pd = build_pipeline(view)
    times = available_times(htg)

    for group_name, time in selected_times(times):
        output_path = args.output_dir / output_filename(
            args.output_prefix,
            group_name,
            time,
        )
        print(f"Rendering {group_name}: t = {time:.12g} -> {output_path}")
        render_time(view, (htg, pd), time, output_path)


if __name__ == "__main__":
    main()
