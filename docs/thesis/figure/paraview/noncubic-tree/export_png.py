"""Render the non-cubic tree ParaView view to PNG.

Run from ParaView's Python shell with an active source, or with pvpython:

    pvpython export_png.py --input case/karman_0.vtkhdf

By default this writes and crops ``../../raster/noncubic-tree_surface_edges.png``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from paraview.simple import (
    ColorBy,
    GetActiveSource,
    GetActiveViewOrCreate,
    GetColorTransferFunction,
    GetLayout,
    GetRepresentation,
    OpenDataFile,
    Render,
    SaveScreenshot,
    Show,
    _DisableFirstRenderCameraReset,
)
from vtkmodules.vtkIOImage import vtkPNGReader, vtkPNGWriter
from vtkmodules.vtkImagingCore import vtkImageClip


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "case" / "karman_0.vtkhdf"
DEFAULT_OUTPUT_DIR = ROOT.parent.parent / "raster"
DEFAULT_OUTPUT_NAME = f"{ROOT.name}_surface_edges.png"

MULT = 8
IMAGE_RESOLUTION = (1620*MULT, 1133*MULT)
VIEW_SIZE = (1620, 1133)
BACKGROUND = (1.0, 1.0, 1.0)

COLOR_FIELD = "f"
EDGE_COLOR = [0.0, 0.0, 0.0]
EDGE_OPACITY = 1.0
EDGE_LINE_WIDTH = 0.25*MULT
CROP_PADDING = 0
ALPHA_THRESHOLD = 0

CAMERA_POSITION = [3.470919377466497, -0.026404926087655968, 26.8]
CAMERA_FOCAL_POINT = [3.470919377466497, -0.026404926087655968, 0.0]
CAMERA_PARALLEL_SCALE = 3.3315114662390672


def parse_args():
    parser = argparse.ArgumentParser(
        description="Render the non-cubic tree surface-with-edges view.",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=(
            "VTKHDF input file to open. If this file does not exist, the "
            "currently active ParaView source is used instead "
            f"(default: {DEFAULT_INPUT})."
        ),
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
        "--output-prefix",
        default="",
        help="Optional filename prefix for script compatibility.",
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
        "--edge-line-width",
        type=float,
        default=EDGE_LINE_WIDTH,
        help=f"Black edge line width (default: {EDGE_LINE_WIDTH}).",
    )
    parser.add_argument(
        "--crop-padding",
        type=int,
        default=CROP_PADDING,
        help=f"Transparent padding to leave after cropping, in pixels (default: {CROP_PADDING}).",
    )
    parser.add_argument(
        "--no-crop",
        action="store_true",
        help="Keep the full screenshot extent instead of cropping transparent margins.",
    )
    return parser.parse_args()


def output_path(args):
    output_name = args.output_name
    if args.output_prefix:
        output_name = f"{args.output_prefix}_{output_name}"
    return args.output_dir / output_name


def get_or_open_source(input_path):
    if input_path.exists():
        return OpenDataFile(str(input_path))

    source = GetActiveSource()
    if source is None:
        raise FileNotFoundError(
            f"No input file found at {input_path} and no active ParaView source exists."
        )
    return source


def configure_view(view_size):
    view = GetActiveViewOrCreate("RenderView")

    layout = GetLayout()
    layout.SetSize(*view_size)

    view.ViewSize = view_size
    view.Background = BACKGROUND
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.Set(
        InteractionMode="2D",
        CameraPosition=CAMERA_POSITION,
        CameraFocalPoint=CAMERA_FOCAL_POINT,
        CameraParallelScale=CAMERA_PARALLEL_SCALE,
    )
    return view


def configure_display(source, view, edge_line_width):
    display = GetRepresentation(source, view=view)
    if display is None:
        display = Show(source, view)

    display.SetRepresentationType("Surface With Edges")
    ColorBy(display, ("CELLS", COLOR_FIELD))
    display.EdgeColor = EDGE_COLOR
    display.EdgeOpacity = EDGE_OPACITY
    display.LineWidth = edge_line_width
    display.SetScalarBarVisibility(view, False)

    GetColorTransferFunction(COLOR_FIELD)
    return display


def save_png(view, path, image_resolution):
    path.parent.mkdir(parents=True, exist_ok=True)
    Render(view)
    SaveScreenshot(
        str(path),
        view,
        ImageResolution=image_resolution,
        TransparentBackground=1,
    )


def crop_png(path, padding):
    reader = vtkPNGReader()
    reader.SetFileName(str(path))
    reader.Update()

    image = reader.GetOutput()
    extent = image.GetExtent()
    scalars = image.GetPointData().GetScalars()
    components = scalars.GetNumberOfComponents()

    if components not in (2, 4):
        print(f"Skipping crop: {path} has no alpha channel.")
        return

    xmin, xmax, ymin, ymax, zmin, _ = extent
    crop_extent = None

    for y in range(ymin, ymax + 1):
        for x in range(xmin, xmax + 1):
            point_id = image.ComputePointId((x, y, zmin))
            alpha = scalars.GetTuple(point_id)[components - 1]
            if alpha <= ALPHA_THRESHOLD:
                continue

            if crop_extent is None:
                crop_extent = [x, x, y, y]
            else:
                crop_extent[0] = min(crop_extent[0], x)
                crop_extent[1] = max(crop_extent[1], x)
                crop_extent[2] = min(crop_extent[2], y)
                crop_extent[3] = max(crop_extent[3], y)

    if crop_extent is None:
        print(f"Skipping crop: {path} is fully transparent.")
        return

    crop_extent[0] = max(xmin, crop_extent[0] - padding)
    crop_extent[1] = min(xmax, crop_extent[1] + padding)
    crop_extent[2] = max(ymin, crop_extent[2] - padding)
    crop_extent[3] = min(ymax, crop_extent[3] + padding)

    clip = vtkImageClip()
    clip.SetInputConnection(reader.GetOutputPort())
    clip.SetOutputWholeExtent(
        crop_extent[0],
        crop_extent[1],
        crop_extent[2],
        crop_extent[3],
        zmin,
        zmin,
    )
    clip.ClipDataOn()
    clip.Update()

    writer = vtkPNGWriter()
    writer.SetFileName(str(path))
    writer.SetInputConnection(clip.GetOutputPort())
    writer.Write()

    width = crop_extent[1] - crop_extent[0] + 1
    height = crop_extent[3] - crop_extent[2] + 1
    print(f"Cropped PNG to {width}x{height} -> {path}")


def main():
    args = parse_args()
    _DisableFirstRenderCameraReset()

    source = get_or_open_source(args.input)
    view = configure_view(args.view_size)
    configure_display(source, view, args.edge_line_width)

    path = output_path(args)
    print(
        f"Rendering ParaView PNG at {args.resolution[0]}x{args.resolution[1]} "
        f"-> {path}"
    )
    save_png(view, path, args.resolution)
    if not args.no_crop:
        crop_png(path, args.crop_padding)


if __name__ == "__main__":
    main()
