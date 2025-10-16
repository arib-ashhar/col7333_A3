import argparse
import sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

DIR_STEP = {
    'L': (-1, 0),
    'R': ( 1, 0),
    'U': ( 0,-1),
    'D': ( 0, 1),
}

def parse_city(path):
    """
    Returns:
      scenario (int),
      N, M, K, J (ints),
      lines: list of (sx,sy, ex,ey),
      popular: list of (x,y) (possibly empty)
    """
    with open(path, 'r') as f:
        tokens = f.read().split()
    it = iter(tokens)

    try:
        scenario = int(next(it))
    except StopIteration:
        raise ValueError("Empty .city file")

    if scenario == 1:
        N = int(next(it)); M = int(next(it)); K = int(next(it)); J = int(next(it))
        lines = []
        for _ in range(K):
            sx = int(next(it)); sy = int(next(it)); ex = int(next(it)); ey = int(next(it))
            lines.append((sx,sy,ex,ey))
        popular = []
    elif scenario == 2:
        N = int(next(it)); M = int(next(it)); K = int(next(it)); J = int(next(it)); P = int(next(it))
        lines = []
        for _ in range(K):
            sx = int(next(it)); sy = int(next(it)); ex = int(next(it)); ey = int(next(it))
            lines.append((sx,sy,ex,ey))
        popular = []
        for _ in range(P):
            px = int(next(it)); py = int(next(it))
            popular.append((px,py))
    else:
        raise ValueError(f"Unsupported scenario type: {scenario}")

    return scenario, N, M, K, J, lines, popular

def parse_metromap(path, K_expected):
    """
    Each line: tokens of directions ending with 0 (e.g., 'R R D 0').
    Returns list[list[str]] of directions per line.
    If file contains a single '0', returns [] and unsat=True.
    """
    dirs = []
    with open(path, 'r') as f:
        raw_lines = [ln.strip() for ln in f if ln.strip()]

    if len(raw_lines) == 1 and raw_lines[0] == '0':
        return [], True

    for ln in raw_lines:
        toks = ln.split()
        if not toks or toks[-1] != '0':
            raise ValueError(f"Bad metromap row (must end with 0): {ln}")
        steps = [t for t in toks[:-1] if t in DIR_STEP]
        dirs.append(steps)

    if K_expected is not None and len(dirs) != K_expected:
        print(f"[warn] metromap has {len(dirs)} rows but K={K_expected}", file=sys.stderr)

    return dirs, False

def build_paths(N, M, lines, dir_rows):
    """
    From starts and direction rows, reconstruct cell sequences.
    Returns list[list[(x,y)]]
    """
    paths = []
    for idx, ((sx,sy,ex,ey), steps) in enumerate(zip(lines, dir_rows)):
        x, y = sx, sy
        seq = [(x,y)]
        for d in steps:
            dx, dy = DIR_STEP[d]
            x, y = x + dx, y + dy
            if x < 0 or x >= N or y < 0 or y >= M:
                raise ValueError(f"Path for line {idx} goes out of bounds at {(x,y)}")
            seq.append((x,y))
        # Optional sanity: end matches
        # If you want a strict check, uncomment:
        # if seq[-1] != (ex,ey):
        #     print(f"[warn] line {idx} ends at {seq[-1]} but expected {(ex,ey)}", file=sys.stderr)
        paths.append(seq)
    return paths

def draw_grid(ax, N, M):
    # draw light grid lines
    ax.set_xlim(-0.5, N-0.5)
    ax.set_ylim(-0.5, M-0.5)
    ax.set_xticks(range(N))
    ax.set_yticks(range(M))
    ax.grid(True, which='both', color='#cccccc', linewidth=0.5)
    ax.set_aspect('equal', adjustable='box')
    ax.invert_yaxis()  # (0,0) at top-left like the spec
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis='both', which='both', labelsize=8)

def plot_map(N, M, K, lines, popular, paths, out_path=None, title="Metro Map Visualization"):
    fig, ax = plt.subplots(figsize=(10,10))
    draw_grid(ax, N, M)
    ax.set_title(title, fontsize=16, pad=14)

    # color cycle per line (repeat if needed)
    colors = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    legend_handles = []
    # Popular cells (stars)
    if popular:
        xs = [x for x,_ in popular]
        ys = [y for _,y in popular]
        ax.scatter(xs, ys, marker='*', s=250, c='gold', edgecolor='k', linewidths=0.8, zorder=5, label='Popular Cells')
        legend_handles.append(Line2D([0],[0], marker='*', color='w', label='Popular Cells',
                                     markerfacecolor='gold', markeredgecolor='k', markersize=15))

    for k, (path, line) in enumerate(zip(paths, lines)):
        col = colors[k % len(colors)]
        # draw polyline
        xs = [x for x,_ in path]
        ys = [y for _,y in path]
        ax.plot(xs, ys, '-', color=col, linewidth=3, alpha=0.9, zorder=2)

        # start/end squares
        sx, sy, ex, ey = line
        ax.scatter([sx], [sy], marker='s', s=120, color=col, edgecolor='k', linewidths=0.8, zorder=6)
        ax.scatter([ex], [ey], marker='s', s=120, color=col, edgecolor='k', linewidths=0.8, zorder=6)

        legend_handles.append(Line2D([0],[0], color=col, lw=3, label=f"Line {k}"))

    # Legend
    if legend_handles:
        ax.legend(handles=legend_handles, loc='center', bbox_to_anchor=(0.5, -0.08),
                  ncol=3, fontsize=9, frameon=True)

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=200)
        print(f"Saved figure to {out_path}")
    plt.show()

def main():
    ap = argparse.ArgumentParser(description="Visualize metro map from .city and .metromap.")
    ap.add_argument("city", help="Input .city file (scenario 1 or 2)")
    ap.add_argument("metromap", help="Input .metromap file (directions)")
    ap.add_argument("-o", "--output", help="Output image path (e.g., out.png)")
    args = ap.parse_args()

    scenario, N, M, K, J, lines, popular = parse_city(args.city)
    dir_rows, unsat = parse_metromap(args.metromap, K_expected=K)
    if unsat:
        print("Instance is UNSAT (metromap contains only 0). Nothing to draw.")
        sys.exit(0)

    # Rebuild paths from starts + directions
    paths = build_paths(N, M, lines, dir_rows)

    title = "Metro Map Visualization"
    if scenario == 2:
        title += " (Scenario 2)"

    # Default output name
    out_path = args.output
    if out_path is None:
        base = Path(args.metromap).with_suffix('')
        out_path = f"{base}_viz.png"

    plot_map(N, M, K, lines, popular, paths, out_path=out_path, title=title)

if __name__ == "__main__":
    main()