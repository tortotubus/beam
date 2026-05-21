"""Compute flapping amplitude and Strouhal number from tip.csv.

The default input is the local ParaView case output:

    python analyze_tip.py

Use ``--csv`` to analyze another file with columns ``t`` and ``y_tip``.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parent
DEFAULT_CSV = ROOT / "case" / "tip.csv"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute tip oscillation amplitude and Strouhal number.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help=f"CSV file with t and y_tip columns (default: {DEFAULT_CSV})",
    )
    parser.add_argument(
        "--start-time",
        type=float,
        default=10.0,
        help="Ignore samples before this nondimensional time.",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=1001,
        help="Centered moving-average window in samples. Use 1 to disable.",
    )
    parser.add_argument(
        "--min-separation",
        type=float,
        default=1.5,
        help="Minimum nondimensional time between extrema of the same type.",
    )
    parser.add_argument(
        "--prominence",
        type=float,
        default=0.05,
        help="Minimum absolute displacement for accepted extrema.",
    )
    parser.add_argument(
        "--keep-first-cycle",
        action="store_true",
        help="Include the first detected cycle after --start-time.",
    )
    return parser.parse_args()


def read_tip_data(csv_path):
    times = []
    values = []

    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required_columns = {"t", "y_tip"}
        missing_columns = required_columns - set(reader.fieldnames or [])
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"{csv_path} is missing required column(s): {missing}")

        for row in reader:
            times.append(float(row["t"]))
            values.append(float(row["y_tip"]))

    return times, values


def moving_average(values, window):
    if window <= 1:
        return list(values)

    if window % 2 == 0:
        raise ValueError("--smooth-window must be odd")

    half_window = window // 2
    prefix = [0.0]
    for value in values:
        prefix.append(prefix[-1] + value)

    smoothed = []
    for index in range(len(values)):
        begin = max(0, index - half_window)
        end = min(len(values), index + half_window + 1)
        smoothed.append((prefix[end] - prefix[begin]) / (end - begin))

    return smoothed


def extrema_candidates(times, values, start_time):
    peaks = []
    troughs = []

    for index in range(1, len(values) - 1):
        if times[index] < start_time:
            continue

        previous_value = values[index - 1]
        value = values[index]
        next_value = values[index + 1]

        if value > previous_value and value >= next_value:
            peaks.append((times[index], value))
        if value < previous_value and value <= next_value:
            troughs.append((times[index], value))

    return peaks, troughs


def select_extrema(candidates, *, mode, min_separation, prominence):
    selected = []

    for time, value in candidates:
        if mode == "peak" and value < prominence:
            continue
        if mode == "trough" and value > -prominence:
            continue

        if selected and time - selected[-1][0] < min_separation:
            previous_time, previous_value = selected[-1]
            if mode == "peak" and value > previous_value:
                selected[-1] = (time, value)
            if mode == "trough" and value < previous_value:
                selected[-1] = (time, value)
        else:
            selected.append((time, value))

    return selected


def periods_from_extrema(peaks, troughs):
    peak_periods = [
        next_peak[0] - peak[0]
        for peak, next_peak in zip(peaks, peaks[1:])
    ]
    trough_periods = [
        next_trough[0] - trough[0]
        for trough, next_trough in zip(troughs, troughs[1:])
    ]
    return peak_periods + trough_periods


def compute_metrics(times, values, args):
    smoothed_values = moving_average(values, args.smooth_window)
    peak_candidates, trough_candidates = extrema_candidates(
        times,
        smoothed_values,
        args.start_time,
    )

    peaks = select_extrema(
        peak_candidates,
        mode="peak",
        min_separation=args.min_separation,
        prominence=args.prominence,
    )
    troughs = select_extrema(
        trough_candidates,
        mode="trough",
        min_separation=args.min_separation,
        prominence=args.prominence,
    )

    if len(peaks) < 2 or len(troughs) < 2:
        raise ValueError("not enough extrema found; adjust start time or thresholds")

    periods = periods_from_extrema(peaks, troughs)
    peak_values = [value for _, value in peaks]
    trough_values = [value for _, value in troughs]

    if not args.keep_first_cycle and len(periods) > 2:
        periods = periods[1:]
        peak_values = peak_values[1:]
        trough_values = trough_values[1:]

    mean_peak = mean(peak_values)
    mean_trough = mean(trough_values)
    amplitude = 0.5 * (mean_peak - mean_trough)
    midline = 0.5 * (mean_peak + mean_trough)
    period = mean(periods)
    strouhal = 1.0 / period

    return {
        "peaks": peaks,
        "troughs": troughs,
        "mean_peak": mean_peak,
        "mean_trough": mean_trough,
        "midline": midline,
        "amplitude": amplitude,
        "period": period,
        "strouhal": strouhal,
    }


def print_metrics(metrics, args):
    print(f"csv: {args.csv}")
    print(f"analysis window: t >= {args.start_time:g}")
    print(f"smoothing window: {args.smooth_window} samples")
    print(f"peaks: {len(metrics['peaks'])}")
    print(f"troughs: {len(metrics['troughs'])}")
    print(f"mean peak y_tip: {metrics['mean_peak']:.6f}")
    print(f"mean trough y_tip: {metrics['mean_trough']:.6f}")
    print(f"midline: {metrics['midline']:.6f}")
    print(f"amplitude: {metrics['amplitude']:.6f}")
    print(f"period: {metrics['period']:.6f}")
    print(f"St: {metrics['strouhal']:.6f}")
    print(
        "peak times: "
        + ", ".join(f"{time:.3f}" for time, _ in metrics["peaks"])
    )
    print(
        "trough times: "
        + ", ".join(f"{time:.3f}" for time, _ in metrics["troughs"])
    )


def main():
    args = parse_args()
    times, values = read_tip_data(args.csv)
    metrics = compute_metrics(times, values, args)
    print_metrics(metrics, args)


if __name__ == "__main__":
    main()
