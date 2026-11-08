#!/usr/bin/env python3
"""
Generates the FBGA wave-crest mark from parametric curve geometry.

Model: a rise phase whose slope falls linearly from RISE_SLOPE_START to
RISE_SLOPE_END over RISE_WIDTH, followed by a brake phase whose slope rises
linearly from BRAKE_SLOPE_START to BRAKE_SLOPE_END over a width solved so the
net change exactly cancels the rise's net gain (same start/end height).
Since y is quadratic in a parameter that's linear in x, each phase maps to a
single SVG quadratic Bezier (`Q` command).

Usage:
    python3 generate_logo.py
    python3 generate_logo.py --rise-start 2 --rise-end 0.1 --brake-start -3 --brake-end -1

Outputs (next to this script):
    fbga-mark.svg            plain icon -- safe down to ~16px
    fbga-badge-carved.svg    same shape with "FBGA" knocked out -- needs ~64px+
    fbga-lockup.svg          icon + wordmark, side by side
    *.png                    rasterized at a few common sizes, if rsvg-convert
                              (or inkscape) is on PATH
"""
import argparse
import math
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

BLUE = "#1E90FF"
ORANGE = "#FF8C00"
DISPLAY_FONT = "Helvetica Neue Condensed, Arial Narrow, Arial, sans-serif"


def quad_bezier_from_slope(width, slope_start, slope_end, v_start):
    """One phase: slope(t) linear in t over `width`, v(0)=v_start.

    Returns (v_end, control_y) where control_y is the SVG Q control point's
    y-value for a Bezier whose x runs linearly 0..width and whose y matches
    the exact quadratic antiderivative of slope(t).
    """
    net = (slope_start + slope_end) / 2 * width
    v_end = v_start + net
    # V(t) = v_start + slope_start*width*t + (slope_end-slope_start)*width/2 * t^2
    a1 = slope_start * width
    a2 = (slope_end - slope_start) * width / 2
    control_y = v_start + a1 / 2
    assert abs((v_start + a1 + a2) - v_end) < 1e-6
    return v_end, control_y


def solve_brake_width(rise_width, rise_net, brake_start, brake_end):
    """Width needed so the brake phase's net change is exactly -rise_net."""
    avg = (brake_start + brake_end) / 2
    if avg == 0:
        raise ValueError("brake average slope is zero; curve can't rejoin the baseline")
    return -rise_net / avg


def build_geometry(rise_start, rise_end, brake_start, brake_end,
                    rise_width=100.0, baseline_ratio=0.088, pad=10.0):
    rise_net = (rise_start + rise_end) / 2 * rise_width
    brake_width = solve_brake_width(rise_width, rise_net, brake_start, brake_end)

    # baseline_ratio = v0 / peak  =>  v0 = ratio * (v0 + rise_net)
    v0 = baseline_ratio * rise_net / (1 - baseline_ratio)
    peak = v0 + rise_net

    _, rise_ctrl_v = quad_bezier_from_slope(rise_width, rise_start, rise_end, v0)
    v_end, brake_ctrl_v = quad_bezier_from_slope(brake_width, brake_start, brake_end, peak)

    # screen space: x = pad + s ; y = (pad + peak) - v   (so v=0 sits on the ground line)
    ground_y = pad + peak
    def P(s, v):
        return (pad + s, ground_y - v)

    start = P(0, v0)
    rise_ctrl = P(rise_width / 2, rise_ctrl_v)
    peak_pt = P(rise_width, peak)
    brake_ctrl = P(rise_width + brake_width / 2, brake_ctrl_v)
    end = P(rise_width + brake_width, v_end)

    width = pad * 2 + rise_width + brake_width
    height = ground_y + pad * 0.3  # small breathing room below the ground line

    return {
        "start": start, "rise_ctrl": rise_ctrl, "peak": peak_pt,
        "brake_ctrl": brake_ctrl, "end": end,
        "width": width, "height": height, "ground_y": ground_y,
        "v0": v0, "peak_v": peak, "rise_width": rise_width, "brake_width": brake_width,
    }


def fmt(n):
    return f"{n:.2f}".rstrip("0").rstrip(".")


def paths(geo):
    s, rc, pk, bc, e, gy = (geo["start"], geo["rise_ctrl"], geo["peak"],
                             geo["brake_ctrl"], geo["end"], geo["ground_y"])
    rise_curve = f"M{fmt(s[0])},{fmt(s[1])} Q{fmt(rc[0])},{fmt(rc[1])} {fmt(pk[0])},{fmt(pk[1])}"
    brake_curve = f"M{fmt(pk[0])},{fmt(pk[1])} Q{fmt(bc[0])},{fmt(bc[1])} {fmt(e[0])},{fmt(e[1])}"
    blue_fill = (f"M{fmt(s[0])},{fmt(gy)} L{fmt(s[0])},{fmt(s[1])} "
                 f"Q{fmt(rc[0])},{fmt(rc[1])} {fmt(pk[0])},{fmt(pk[1])} "
                 f"L{fmt(pk[0])},{fmt(gy)} Z")
    orange_fill = (f"M{fmt(pk[0])},{fmt(gy)} L{fmt(pk[0])},{fmt(pk[1])} "
                   f"Q{fmt(bc[0])},{fmt(bc[1])} {fmt(e[0])},{fmt(e[1])} "
                   f"L{fmt(e[0])},{fmt(gy)} Z")
    return rise_curve, brake_curve, blue_fill, orange_fill


def write_mark_svg(geo, out_path):
    _, _, blue_fill, orange_fill = paths(geo)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {fmt(geo['width'])} {fmt(geo['height'])}">
  <path fill="{BLUE}" d="{blue_fill}"/>
  <path fill="{ORANGE}" d="{orange_fill}"/>
</svg>
'''
    out_path.write_text(svg)


def write_carved_badge_svg(geo, out_path, letters=("F", "B", "G", "A"),
                            font_size=36, baseline_y=None, letter_gap=30):
    _, _, blue_fill, orange_fill = paths(geo)
    w, h = geo["width"], geo["height"]
    if baseline_y is None:
        baseline_y = geo["ground_y"] - 6  # a hair above the true ground

    # one uniform gap across all four letters, centered on the whole mark --
    # NOT per-color pairs, which read as two separate words ("FB" "GA").
    mark_center = (geo["start"][0] + geo["end"][0]) / 2
    centers = [mark_center + k * letter_gap for k in (-1.5, -0.5, 0.5, 1.5)]

    text_nodes = "\n    ".join(
        f'<text x="{fmt(cx)}" y="{fmt(baseline_y)}" text-anchor="middle" '
        f'font-family="{DISPLAY_FONT}" font-weight="800" font-size="{font_size}" fill="black">{ch}</text>'
        for cx, ch in zip(centers, letters)
    )

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {fmt(w)} {fmt(h)}">
  <defs>
    <mask id="fbga-cutout" maskUnits="userSpaceOnUse" x="0" y="0" width="{fmt(w)}" height="{fmt(h)}">
      <rect x="0" y="0" width="{fmt(w)}" height="{fmt(h)}" fill="white"/>
      {text_nodes}
    </mask>
  </defs>
  <g mask="url(#fbga-cutout)">
    <path fill="{BLUE}" d="{blue_fill}"/>
    <path fill="{ORANGE}" d="{orange_fill}"/>
  </g>
</svg>
'''
    out_path.write_text(svg)


def write_lockup_svg(geo, out_path):
    _, _, blue_fill, orange_fill = paths(geo)
    w, h = geo["width"], geo["height"]
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {fmt(w*2.4)} {fmt(h)}">
  <g transform="translate(0,0)">
    <path fill="{BLUE}" d="{blue_fill}"/>
    <path fill="{ORANGE}" d="{orange_fill}"/>
  </g>
  <text x="{fmt(w*1.15)}" y="{fmt(h*0.78)}" font-family="{DISPLAY_FONT}" font-weight="800"
        font-size="{fmt(h*0.62)}" letter-spacing="2">FBGA</text>
</svg>
'''
    out_path.write_text(svg)


def export_pngs(svg_path, sizes):
    """Export at each `size` as the longer side, aspect ratio preserved.

    -a keeps the real proportions instead of forcing a square canvas, which
    would stretch anything that isn't already square.
    """
    tool = shutil.which("rsvg-convert") or shutil.which("inkscape")
    if not tool:
        print(f"  (skipping PNG export for {svg_path.name} -- no rsvg-convert/inkscape on PATH)")
        return
    for size in sizes:
        out = svg_path.with_name(f"{svg_path.stem}-{size}.png")
        if "rsvg-convert" in tool:
            cmd = [tool, "-w", str(size), "-h", str(size), "-a", str(svg_path), "-o", str(out)]
        else:
            cmd = [tool, str(svg_path), f"--export-filename={out}", f"-w={size}"]
        subprocess.run(cmd, check=True)
        print(f"  wrote {out.name}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rise-start", type=float, default=2.0)
    ap.add_argument("--rise-end", type=float, default=0.1)
    ap.add_argument("--brake-start", type=float, default=-2.0)
    ap.add_argument("--brake-end", type=float, default=-0.1)
    ap.add_argument("--baseline-ratio", type=float, default=0.088,
                     help="start/end velocity as a fraction of the peak (0.05-0.10 looks right)")
    ap.add_argument("--rise-width", type=float, default=100.0)
    ap.add_argument("--no-png", action="store_true", help="skip PNG export even if a rasterizer is available")
    args = ap.parse_args()

    geo = build_geometry(args.rise_start, args.rise_end, args.brake_start, args.brake_end,
                          rise_width=args.rise_width, baseline_ratio=args.baseline_ratio)

    print(f"rise width={geo['rise_width']:.1f}  brake width={geo['brake_width']:.1f}  "
          f"(brake/rise = {geo['brake_width']/geo['rise_width']*100:.1f}%)")
    print(f"v0={geo['v0']:.2f}  peak={geo['peak_v']:.2f}  baseline ratio={geo['v0']/geo['peak_v']*100:.1f}%")

    mark_path = HERE / "fbga-mark.svg"
    badge_path = HERE / "fbga-badge-carved.svg"
    lockup_path = HERE / "fbga-lockup.svg"

    write_mark_svg(geo, mark_path)
    write_carved_badge_svg(geo, badge_path)
    write_lockup_svg(geo, lockup_path)
    print(f"wrote {mark_path.name}, {badge_path.name}, {lockup_path.name}")

    if not args.no_png:
        export_pngs(mark_path, [16, 32, 64, 192, 512])
        export_pngs(badge_path, [64, 128, 256, 512])


if __name__ == "__main__":
    main()
