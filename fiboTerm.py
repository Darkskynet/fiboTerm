#!/usr/bin/env python3
"""
fiboTerm.py

A terminal "deluxe" golden spiral visualiser.

Features:
 - Single continuous golden/logarithmic spiral arm (radius grows by φ every 90°)
 - Auto-detect terminal size and scale spiral to fit
 - --braille: higher resolution rendering using Unicode Braille (2x4 subpixels per char)
 - --trail: fading trail
 - --spin: rotate the completed spiral in-place
 - --pulse: rhythmic brightness pulse
 - --color [rainbow|pulse|none] with --color-source [time|distance]
 - --goldenratio: display φ
 - --fibs: print fibonacci numbers used
 - --hypnotic: repeat drawing forever
 - --turns, --delay to control drawing
 - No external dependencies (stdlib only)

Example:
    python fiboTerm.py --turns 8 --braille --trail --color rainbow --color-source distance
"""

import argparse
import math
import os
import sys
import time
import shutil
from typing import Tuple, List, Optional

# ---------------------------------------------------------------------------
# Utility & constants
# ---------------------------------------------------------------------------
PHI = (1 + 5 ** 0.5) / 2.0
CSI = "\x1b["
RESET = "\x1b[0m"

# Braille helper constants
BRAILLE_BASE = 0x2800
# Dot bit mapping in braille cell:
# dots:
# 1 4
# 2 5
# 3 6
# 7 8
# bits indices (0-based) correspond to dot numbers - 1 except dot7->6, dot8->7
# mapping: (col, row) -> bit index
BRAILLE_BIT_MAP = {
    (0, 0): 0,  # dot1
    (0, 1): 1,  # dot2
    (0, 2): 2,  # dot3
    (0, 3): 6,  # dot7
    (1, 0): 3,  # dot4
    (1, 1): 4,  # dot5
    (1, 2): 5,  # dot6
    (1, 3): 7,  # dot8
}

# ---------------------------------------------------------------------------
# Terminal helpers
# ---------------------------------------------------------------------------
def get_term_size() -> Tuple[int, int]:
    sz = shutil.get_terminal_size(fallback=(80, 40))
    return sz.columns, sz.lines

def hide_cursor():
    sys.stdout.write(CSI + "?25l")
    sys.stdout.flush()

def show_cursor():
    sys.stdout.write(CSI + "?25h")
    sys.stdout.flush()

def clear_screen():
    sys.stdout.write(CSI + "2J")
    sys.stdout.write(CSI + "H")
    sys.stdout.flush()

# ---------------------------------------------------------------------------
# Color helpers (convert HSV -> xterm-256 color index)
# ---------------------------------------------------------------------------
def hsv_to_rgb(h: float, s: float, v: float) -> Tuple[float, float, float]:
    """Return r,g,b in [0..1]. h in degrees 0..360"""
    h = h % 360
    c = v * s
    x = c * (1 - abs((h / 60.0) % 2 - 1))
    m = v - c
    if h < 60:
        r, g, b = c, x, 0
    elif h < 120:
        r, g, b = x, c, 0
    elif h < 180:
        r, g, b = 0, c, x
    elif h < 240:
        r, g, b = 0, x, c
    elif h < 300:
        r, g, b = x, 0, c
    else:
        r, g, b = c, 0, x
    return r + m, g + m, b + m

def rgb_to_xterm256(r: float, g: float, b: float) -> int:
    """
    Map r,g,b in [0,1] to a 256-color xterm index (16-231 color cube) approximated.
    This is good enough for smooth rainbow transitions in terminals that support 256 colors.
    """
    # Scale to 0..5
    rc = int(round(r * 5))
    gc = int(round(g * 5))
    bc = int(round(b * 5))
    return 16 + 36 * rc + 6 * gc + bc

def hue_color_index(hue_deg: float, sat: float = 1.0, val: float = 1.0) -> int:
    r, g, b = hsv_to_rgb(hue_deg, sat, val)
    return rgb_to_xterm256(r, g, b)

def ansi_fg_xterm(idx: int) -> str:
    return f"\x1b[38;5;{idx}m"

# ---------------------------------------------------------------------------
# Braille grid and ASCII grid classes
# ---------------------------------------------------------------------------
class GridBase:
    def __init__(self, cols: int, rows: int):
        self.cols = cols
        self.rows = rows
        # each cell stores an intensity/fade value (float)
        self.buf = [[0.0 for _ in range(cols)] for _ in range(rows)]
        # char buffer (for ASCII mode) or braille dotmask buffer (for braille)
        self.charbuf = [[" " for _ in range(cols)] for _ in range(rows)]

    def clear(self):
        for y in range(self.rows):
            for x in range(self.cols):
                self.buf[y][x] = 0.0
                self.charbuf[y][x] = " "

    def decay(self, amount: float):
        for y in range(self.rows):
            for x in range(self.cols):
                self.buf[y][x] = max(0.0, self.buf[y][x] - amount)
                if self.buf[y][x] == 0.0:
                    self.charbuf[y][x] = " "

class AsciiGrid(GridBase):
    def __init__(self, cols: int, rows: int):
        super().__init__(cols, rows)

    def set_pixel(self, x: int, y: int, intensity: float, ch: str):
        if 0 <= x < self.cols and 0 <= y < self.rows:
            if intensity >= self.buf[y][x]:
                self.buf[y][x] = intensity
                self.charbuf[y][x] = ch

    def render_to_screen(self, color_map_func):
        out_lines = []
        for y in range(self.rows):
            line = []
            for x in range(self.cols):
                intensity = self.buf[y][x]
                ch = self.charbuf[y][x] if intensity > 0 else " "
                if intensity > 0:
                    color_code = color_map_func(x, y, intensity)
                    line.append(f"{ansi_fg_xterm(color_code)}{ch}{RESET}")
                else:
                    line.append(" ")
            out_lines.append("".join(line))
        sys.stdout.write(CSI + "H")
        sys.stdout.write("\n".join(out_lines))
        sys.stdout.flush()

class BrailleGrid:
    """
    A higher resolution grid using braille cells. Each braille char cell is 2x4 subpixels.
    We store a high-resolution boolean/intensity canvas and then compress to braille for rendering.
    """
    def __init__(self, char_cols: int, char_rows: int):
        self.char_cols = char_cols
        self.char_rows = char_rows
        self.width = char_cols * 2     # subpixels width
        self.height = char_rows * 4    # subpixels height
        # store intensity per subpixel
        self.sub = [[0.0 for _ in range(self.width)] for _ in range(self.height)]
        # cache braille char mask per cell
        self.charbuf = [[" " for _ in range(self.char_cols)] for _ in range(self.char_rows)]

    def clear(self):
        for y in range(self.height):
            for x in range(self.width):
                self.sub[y][x] = 0.0
        for y in range(self.char_rows):
            for x in range(self.char_cols):
                self.charbuf[y][x] = " "

    def decay(self, amount: float):
        for y in range(self.height):
            for x in range(self.width):
                self.sub[y][x] = max(0.0, self.sub[y][x] - amount)

    def set_subpixel(self, sx: int, sy: int, intensity: float):
        if 0 <= sx < self.width and 0 <= sy < self.height:
            self.sub[sy][sx] = max(self.sub[sy][sx], intensity)

    def compress_to_braille(self):
        # produce charbuf with braille unicode characters
        for cy in range(self.char_rows):
            for cx in range(self.char_cols):
                mask = 0
                # iterate over subpixel block 2x4
                for sub_y in range(4):
                    for sub_x in range(2):
                        sx = cx * 2 + sub_x
                        sy = cy * 4 + sub_y
                        if self.sub[sy][sx] > 0.0:
                            bit = BRAILLE_BIT_MAP[(sub_x, sub_y)]
                            mask |= 1 << bit
                if mask == 0:
                    self.charbuf[cy][cx] = " "
                else:
                    ch = chr(BRAILLE_BASE + mask)
                    self.charbuf[cy][cx] = ch

    def render_to_screen(self, color_map_func):
        self.compress_to_braille()
        out_lines = []
        for y in range(self.char_rows):
            line = []
            for x in range(self.char_cols):
                ch = self.charbuf[y][x]
                if ch != " ":
                    # Determine mean intensity of the block for coloring
                    total = 0.0
                    count = 0
                    for sub_y in range(4):
                        for sub_x in range(2):
                            sx = x * 2 + sub_x
                            sy = y * 4 + sub_y
                            total += self.sub[sy][sx]
                            count += 1
                    mean_int = total / count if count else 0.0
                    color_code = color_map_func(x, y, mean_int)
                    line.append(f"{ansi_fg_xterm(color_code)}{ch}{RESET}")
                else:
                    line.append(" ")
            out_lines.append("".join(line))
        sys.stdout.write(CSI + "H")
        sys.stdout.write("\n".join(out_lines))
        sys.stdout.flush()

# ---------------------------------------------------------------------------
# Spiral generation (single continuous golden/logarithmic spiral)
# ---------------------------------------------------------------------------
def generate_golden_spiral_points(turns: int, points_per_degree: int = 4) -> List[Tuple[float, float, float]]:
    """
    Generate points (x,y,theta_degrees) along a logarithmic golden spiral such that
    r(theta) = a * phi^(theta/90°). We return a list of (x,y,deg) up to 90*turns degrees.
    """
    a = 1.0
    total_degrees = 90 * turns
    steps = max(2, int(total_degrees * points_per_degree))
    pts = []
    for s in range(steps + 1):
        deg = s * (total_degrees / steps)
        theta = math.radians(deg)
        r = a * (PHI ** (deg / 90.0))
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        pts.append((x, y, deg))
    return pts

# ---------------------------------------------------------------------------
# ASCII direction chooser
# ---------------------------------------------------------------------------
def choose_ascii_char(dx: float, dy: float) -> str:
    angle = math.degrees(math.atan2(dy, dx)) % 180
    # map ranges to characters (makes curves look smoother)
    if 67.5 <= angle < 112.5:
        return "|"
    elif 22.5 <= angle < 67.5:
        return "/"
    elif 112.5 <= angle < 157.5:
        return "\\"
    else:
        return "-"

# ---------------------------------------------------------------------------
# Color mapping wrapper: returns xterm256 index based on mode and source
# ---------------------------------------------------------------------------
class Colorizer:
    def __init__(self, mode: str = "rainbow", source: str = "distance", pulse: bool = False):
        self.mode = mode  # rainbow | pulse | none
        self.source = source  # time | distance
        self.pulse = pulse
        self.t0 = time.time()

    def color_for(self, x: int, y: int, intensity: float, *,
                  distance_val: float = 0.0) -> int:
        """
        Determine color index for a given grid cell.
        distance_val is a normalized 0..1 value (0 center, 1 far).
        intensity is 0..1 to allow dimming for trail.
        """
        # base hue source
        src = 0.0
        if self.source == "time":
            src = (time.time() - self.t0) % 10.0 / 10.0  # 0..1 repeating every 10s
        else:
            src = distance_val % 1.0

        if self.mode == "none":
            hue = 0.0
            sat = 0.0
            val = max(0.2, min(1.0, intensity))
        elif self.mode == "pulse":
            hue = 200.0  # bluish default
            sat = 0.9
            # pulsate value
            cyc = math.sin((time.time() - self.t0) * 2.0) * 0.5 + 0.5
            val = 0.2 + 0.8 * cyc * intensity
        else:  # rainbow
            hue = src * 360.0
            sat = 0.9
            val = max(0.2, min(1.0, intensity))

            if self.pulse:
                pul = math.sin((time.time() - self.t0) * 2.5) * 0.5 + 0.5
                val = max(0.2, val * (0.5 + 0.5 * pul))

        # convert to color index
        idx = hue_color_index(hue, sat, val)
        return idx

# ---------------------------------------------------------------------------
# Rendering & main loop
# ---------------------------------------------------------------------------
def run(args):
    # get terminal size and compute grid dimensions
    term_w, term_h = get_term_size()
    # keep at least 4 lines for status/instructions
    available_h = max(8, term_h - 4)
    available_w = max(20, term_w)

    # choose braille or ascii grid
    if args.braille:
        # each character is 2x4 subpixels. Reserve full width/height for characters.
        char_cols = available_w
        char_rows = available_h
        grid = BrailleGrid(char_cols, char_rows)
        high_w, high_h = grid.width, grid.height
    else:
        char_cols = available_w
        char_rows = available_h
        grid = AsciiGrid(char_cols, char_rows)
        high_w, high_h = char_cols, char_rows

    # generate spiral points (in mathematical coordinates)
    points = generate_golden_spiral_points(args.turns, points_per_degree=args.density)
    # compute bounds and center them on the grid
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    # scale factor to map spiral coords to subpixel grid.
    # leave a small margin (10%)
    pad = 0.9
    span_x = max_x - min_x if max_x > min_x else 1.0
    span_y = max_y - min_y if max_y > min_y else 1.0

    scale_x = (high_w - 2) * pad / span_x
    scale_y = (high_h - 2) * pad / span_y
    scale = min(scale_x, scale_y)

    # center offset in subpixel coords
    center_x = high_w / 2.0
    center_y = high_h / 2.0

    # precompute mapped subpixel coordinates for every spiral point
    mapped = []
    max_dist = 0.0
    for (x, y, deg) in points:
        sx = (x - (min_x + max_x) / 2.0) * scale + center_x
        sy = (y - (min_y + max_y) / 2.0) * scale + center_y
        mapped.append((sx, sy, deg))
        # distance from center (for coloring)
        d = math.hypot(sx - center_x, sy - center_y)
        if d > max_dist:
            max_dist = d

    # prepare colorizer
    colorizer = Colorizer(mode=args.color, source=args.color_source, pulse=args.pulse)

    # hide cursor and clear
    hide_cursor()
    clear_screen()

    try:
        iteration = 0
        while True:  # supports hypnotic mode repeating
            iteration += 1
            grid.clear()
            tstart = time.time()

            # draw the spiral point-by-point to show growth
            last_sx, last_sy = None, None
            for i, (sx, sy, deg) in enumerate(mapped):
                # normalized distance for color mapping
                dist = math.hypot(sx - center_x, sy - center_y)
                normalized_dist = dist / (max_dist if max_dist > 0 else 1.0)

                if args.braille:
                    # subpixel coords are integer positions on sub grid
                    ix = int(round(sx))
                    iy = int(round(sy))
                    # set a small stroke by painting a few neighboring subpixels
                    intensity = 1.0
                    for dx in (-1, 0, 1):
                        for dy in (-1, 0, 1):
                            grid.set_subpixel(ix + dx, iy + dy, intensity)
                else:
                    # ascii mode: map to char grid
                    gx = int(round(sx))
                    gy = int(round(sy))
                    # pick char by local slope using next point if available
                    if i + 1 < len(mapped):
                        nx, ny, _ = mapped[i + 1]
                        ch = choose_ascii_char(nx - sx, ny - sy)
                    else:
                        ch = "*"
                    grid.set_pixel(gx, gy, 1.0, ch)

                # render incremental frame (or batch render every N points to reduce overhead)
                if (i % args.batch == 0) or (i == len(mapped) - 1):
                    # map color function closure for current painting
                    def color_map(cx, cy, intensity_local, dist_val=normalized_dist):
                        # For braille, cx,cy are character indices; for ascii they're char indices
                        return colorizer.color_for(cx, cy, intensity_local, distance_val=dist_val)

                    grid.render_to_screen(color_map)

                    # status lines below the drawing area
                    sys.stdout.write(CSI + f"{char_rows+1};1H")
                    sys.stdout.write(" " * (char_cols))
                    sys.stdout.write(CSI + f"{char_rows+1};1H")
                    # progress and fibs
                    progress = (i + 1) / len(mapped)
                    barw = min(40, char_cols - 10)
                    filled = int(progress * barw)
                    progbar = "[" + "#" * filled + " " * (barw - filled) + "]"
                    sys.stdout.write(f"Turn {min(args.turns, args.turns)} | {progbar} ")

                    if args.fibs:
                        # show an approximation of Fibonacci values for each 90° step
                        fibs_display = []
                        # compute fibonacci numbers approximated by phi^n / sqrt(5) maybe; but we can display 0..n
                        # Simpler: generate classic fibs up to turns+1
                        fibs = [0, 1]
                        for _ in range(2, args.turns + 5):
                            fibs.append(fibs[-1] + fibs[-2])
                        sys.stdout.write(" Fibs: " + ",".join(str(f) for f in fibs[:args.turns+1]))
                    sys.stdout.flush()

                # optional tiny delay to animate
                if args.delay > 0:
                    time.sleep(args.delay)

                # apply trail decay incrementally
                if args.trail:
                    grid.decay(args.trail_decay_step)
                else:
                    # if not trail, we still keep intensity at max (no decay)
                    pass

            # final render for the completed spiral
            def final_color_map(cx, cy, intensity_local, distance_val=0.0):
                # estimate distance_val by mapping char cell position to distance from center
                if args.braille:
                    # cx,cy are char indices: compute center of the char in subpixel coords
                    cx_sub = cx * 2 + 1
                    cy_sub = cy * 4 + 2
                    dist = math.hypot(cx_sub - center_x, cy_sub - center_y)
                    dv = dist / (max_dist if max_dist > 0 else 1.0)
                else:
                    dist = math.hypot(cx - center_x, cy - center_y)
                    dv = dist / (max_dist if max_dist > 0 else 1.0)
                return colorizer.color_for(cx, cy, intensity_local, distance_val=dv)

            grid.render_to_screen(final_color_map)
            # print golden ratio optionally and other status
            sys.stdout.write(CSI + f"{char_rows+1};1H")
            sys.stdout.write(" " * (char_cols))
            sys.stdout.write(CSI + f"{char_rows+1};1H")
            sys.stdout.write(f"Turns {args.turns} drawn. ")
            if args.goldenratio:
                sys.stdout.write(f"φ ≈ {PHI:.6f}  ")
            sys.stdout.write(f"Mode: color={args.color}, source={args.color_source}")
            sys.stdout.flush()

            # spin effect
            if args.spin:
                do_spin(grid, final_color_map, args.spin_steps, args.spin_delay, center_x, center_y, args.braille)

            # if hypnotic, repeat; else break
            if args.hypnotic:
                # small pause then regenerate with slight parameter adjustments (or same)
                time.sleep(0.5)
                continue
            else:
                break

    finally:
        show_cursor()
        sys.stdout.write(RESET)
        sys.stdout.write("\nDone!\n")
        sys.stdout.flush()

# ---------------------------------------------------------------------------
# Spin implementation: rotate the existing pixel buffer in place and re-render
# (we approximate by sampling the current grid into a rotated target)
# ---------------------------------------------------------------------------
def do_spin(grid_obj, color_map_func, steps: int, delay: float, center_x: float, center_y: float, is_braille: bool):
    # sample buffer into a float intensity map then rotate
    if is_braille:
        w, h = grid_obj.char_cols, grid_obj.char_rows
        sub_w, sub_h = grid_obj.width, grid_obj.height
        # build intensity map from subpixels
        base_map = [[0.0 for _ in range(sub_w)] for _ in range(sub_h)]
        for y in range(sub_h):
            for x in range(sub_w):
                base_map[y][x] = grid_obj.sub[y][x]
        # rotation center in subpixel coords
        cx = center_x
        cy = center_y
    else:
        w, h = grid_obj.cols, grid_obj.rows
        base_map = [[grid_obj.buf[y][x] for x in range(w)] for y in range(h)]
        cx = center_x
        cy = center_y

    for step in range(steps):
        angle = (step / steps) * 2.0 * math.pi
        sin_a = math.sin(angle)
        cos_a = math.cos(angle)

        # create new target buffer
        if is_braille:
            target = [[0.0 for _ in range(sub_w)] for _ in range(sub_h)]
            for ty in range(sub_h):
                for tx in range(sub_w):
                    # map target pixel to source via inverse rotation
                    dx = tx - cx
                    dy = ty - cy
                    sx =  cos_a * dx + sin_a * dy + cx
                    sy = -sin_a * dx + cos_a * dy + cy
                    # sample nearest
                    sx_i = int(round(sx))
                    sy_i = int(round(sy))
                    if 0 <= sx_i < sub_w and 0 <= sy_i < sub_h:
                        target[ty][tx] = base_map[sy_i][sx_i]
            # copy back
            grid_obj.sub = target
            grid_obj.render_to_screen(color_map_func)
        else:
            target = [[0.0 for _ in range(w)] for _ in range(h)]
            for ty in range(h):
                for tx in range(w):
                    dx = tx - cx
                    dy = ty - cy
                    sx =  cos_a * dx + sin_a * dy + cx
                    sy = -sin_a * dx + cos_a * dy + cy
                    sx_i = int(round(sx))
                    sy_i = int(round(sy))
                    if 0 <= sx_i < w and 0 <= sy_i < h:
                        target[ty][tx] = base_map[sy_i][sx_i]
            grid_obj.buf = target
            # derive charbuf from existing char shapes roughly (we keep same char for nonzero)
            for y in range(h):
                for x in range(w):
                    if target[y][x] > 0:
                        # keep previous char if possible, else choose '-'
                        try:
                            ch = grid_obj.charbuf[y][x]
                            if ch == " ":
                                ch = "-"
                        except Exception:
                            ch = "-"
                        grid_obj.charbuf[y][x] = ch
            grid_obj.render_to_screen(color_map_func)

        time.sleep(delay)

# ---------------------------------------------------------------------------
# Arg parsing & entry
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Fibonacci / Golden spiral deluxe visualiser")
    p.add_argument("--braille", action="store_true", help="Use Unicode braille (2x4 subpixels per char)")
    p.add_argument("--trail", action="store_true", help="Enable fading trail")
    p.add_argument("--trail-decay", type=float, default=0.02, dest="trail_decay_step",
                   help="Trail decay per batch render (default 0.02)")
    p.add_argument("--spin", action="store_true", help="Spin the completed spiral in place")
    p.add_argument("--spin-steps", type=int, default=80, help="Spin steps (default 80)")
    p.add_argument("--spin-delay", type=float, default=0.02, dest="spin_delay",
                   help="Delay between spin frames (default 0.02s)")
    p.add_argument("--pulse", action="store_true", help="Enable pulse brightness modulation")
    p.add_argument("--color", choices=["rainbow", "pulse", "none"], default="rainbow", help="Color mode")
    p.add_argument("--color-source", choices=["time", "distance"], default="distance",
                   help="Source for color variation")
    p.add_argument("--goldenratio", action="store_true", help="Display phi value on screen")
    p.add_argument("--fibs", action="store_true", help="Display Fibonacci numbers used")
    p.add_argument("--hypnotic", action="store_true", help="Repeat drawing forever")
    p.add_argument("--turns", type=int, default=6, help="Number of 90° turns to draw (default 6)")
    p.add_argument("--delay", type=float, default=0.0, help="Delay between plotting points (seconds)")
    p.add_argument("--density", type=int, default=4, help="Points per degree density (default 4)")
    p.add_argument("--batch", type=int, default=4, help="Render every N points to improve speed (default 4)")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    # sanity
    try:
        run(args)
    except KeyboardInterrupt:
        show_cursor()
        sys.stdout.write(RESET)
        sys.stdout.write("\nInterrupted.\n")
        sys.exit(0)

