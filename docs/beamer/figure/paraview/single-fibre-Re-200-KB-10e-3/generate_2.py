"""Render the single-fibre Re=200 HTG surface/PD points view at 49 s.

Run with:

    pvpython generate_2.py

By default this writes one PNG into ../../raster.
Use --video to write a numbered PNG series and encode it with ffmpeg.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

from paraview.simple import (
    ColorBy,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
    GetLayout,
    GetOpacityTransferFunction,
    GetScalarBar,
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
DEFAULT_SERIES_PREFIX = f"{ROOT.name}_surface_edges"
DEFAULT_VIDEO_NAME = f"{ROOT.name}_surface_edges.mp4"

HTG_SERIES = CASE_DIR / "htg.vtkhdf.series"
PD_SERIES = CASE_DIR / "pd.vtkhdf.series"

TARGET_TIME = 49.0
IMAGE_RESOLUTION = (2400, 1200)
VIEW_SIZE = (2400, 1200)
BACKGROUND = (1.0, 1.0, 1.0)

HTG_U_FIELD = "u"
HTG_U_COMPONENT = "Magnitude"
HTG_U_RANGE = (0.0, 1.5)
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
    parser.add_argument(
        "--series",
        action="store_true",
        help="Render a numbered PNG frame series instead of one image.",
    )
    parser.add_argument(
        "--series-prefix",
        default=DEFAULT_SERIES_PREFIX,
        help=(
            "Filename prefix for --series frames "
            f"(default: {DEFAULT_SERIES_PREFIX})."
        ),
    )
    parser.add_argument(
        "--series-start-time",
        type=float,
        default=None,
        help="First requested time for --series; default is the first stored time.",
    )
    parser.add_argument(
        "--series-end-time",
        type=float,
        default=None,
        help="Last requested time for --series; default is --time.",
    )
    parser.add_argument(
        "--series-frame-count",
        type=int,
        default=1000,
        help="Number of requested frames for --series before snapping to stored times.",
    )
    parser.add_argument(
        "--video",
        action="store_true",
        help="Render a numbered PNG frame series and encode it to MP4 with ffmpeg.",
    )
    parser.add_argument(
        "--video-output",
        type=Path,
        default=None,
        help=(
            "Output MP4 path for --video "
            f"(default: OUTPUT_DIR/{DEFAULT_VIDEO_NAME})."
        ),
    )
    parser.add_argument(
        "--video-framerate",
        type=float,
        default=24.0,
        help="Frame rate passed to ffmpeg for --video (default: 24).",
    )
    parser.add_argument(
        "--video-crf",
        type=int,
        default=18,
        help="libx264 CRF quality for --video; lower is higher quality (default: 18).",
    )
    parser.add_argument(
        "--ffmpeg-executable",
        default="ffmpeg",
        help="ffmpeg executable name or path used by --video (default: ffmpeg).",
    )
    return parser.parse_args()


def create_view(view_size):
    view = GetActiveViewOrCreate("RenderView")

    layout = GetLayout()
    layout.SetSize(*view_size)

    view.ViewSize = view_size
    view.Background = BACKGROUND
    view.OrientationAxesVisibility = 1
    view.UseColorPaletteForBackground = 0
    view.Set(
        InteractionMode="2D",
        CameraPosition=[0.4, 0, 43.],
        CameraFocalPoint=[0.4, 0, 0.0],
        CameraViewUp=[-1.0, 0.0, 0.0],
        CameraParallelScale=0.75,
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
    htg_display.SetScalarBarVisibility(view, True)

    u_lut = GetColorTransferFunction(HTG_U_FIELD)
    u_lut.RescaleTransferFunction(*HTG_U_RANGE)
    if hasattr(u_lut, "AutomaticRescaleRangeMode"):
        u_lut.AutomaticRescaleRangeMode = "Never"

    u_pwf = GetOpacityTransferFunction(HTG_U_FIELD)
    u_pwf.RescaleTransferFunction(*HTG_U_RANGE)

    scalar_bar = GetScalarBar(u_lut, view)
    scalar_bar.TitleFontSize = 40
    scalar_bar.LabelFontSize = 40
    scalar_bar.ScalarBarThickness = 24
    scalar_bar.WindowLocation = "Any Location"
    scalar_bar.Position = [0.03466204506065853, 0.22981366459627334]
    scalar_bar.ScalarBarLength = 0.7312422360248443

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


def linspace(start, stop, count):
    if count <= 1:
        return [stop]

    step = (stop - start) / (count - 1)
    return [start + index * step for index in range(count)]


def unique_in_order(values):
    unique = []
    for value in values:
        if value not in unique:
            unique.append(value)
    return unique


def selected_series_times(times, args):
    start = (
        nearest_time(times, args.series_start_time)
        if args.series_start_time is not None
        else times[0]
    )
    end = nearest_time(
        times,
        args.series_end_time if args.series_end_time is not None else args.time,
    )

    if end < start:
        start, end = end, start

    requested_times = linspace(start, end, args.series_frame_count)
    return unique_in_order(nearest_time(times, time) for time in requested_times)


def series_output_name(prefix, index):
    return f"{prefix}_{index:04d}.png"


def resolve_executable(executable):
    path = shutil.which(executable)
    if path is not None:
        return path

    explicit_path = Path(executable)
    if explicit_path.exists():
        return str(explicit_path)

    raise FileNotFoundError(
        f"Could not find '{executable}'. Install ffmpeg, add it to PATH, "
        "or pass --ffmpeg-executable /path/to/ffmpeg."
    )


def encode_video(args, ffmpeg):
    output_path = args.video_output
    if output_path is None:
        output_path = args.output_dir / DEFAULT_VIDEO_NAME
    output_path.parent.mkdir(parents=True, exist_ok=True)

    frame_pattern = args.output_dir / f"{args.series_prefix}_%04d.png"
    command = [
        ffmpeg,
        "-y",
        "-framerate",
        f"{args.video_framerate:g}",
        "-i",
        str(frame_pattern),
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        str(args.video_crf),
        str(output_path),
    ]

    print("Encoding video:", " ".join(command))
    subprocess.run(command, check=True)
    return output_path


def render_time(view, sources, htg_display, time, output_path, image_resolution):
    for source in sources:
        source.UpdatePipeline(time)

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

    ffmpeg = None
    if args.video:
        ffmpeg = resolve_executable(args.ffmpeg_executable)

    view = create_view(args.view_size)
    htg, pd, htg_display = build_pipeline(view)
    times = available_times(htg)

    if args.series or args.video:
        selected_times = selected_series_times(times, args)
        for index, time in enumerate(selected_times):
            output_path = args.output_dir / series_output_name(
                args.series_prefix,
                index,
            )
            print(
                f"Rendering frame {index:04d}: t = {time:.12g} "
                f"at {args.resolution[0]}x{args.resolution[1]} -> {output_path}"
            )
            render_time(
                view,
                (htg, pd),
                htg_display,
                time,
                output_path,
                args.resolution,
            )
        if args.video:
            video_path = encode_video(args, ffmpeg)
            print(f"Wrote video: {video_path}")
        return

    time = nearest_time(times, args.time)

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
