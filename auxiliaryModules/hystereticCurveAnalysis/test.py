import numpy as np
import matplotlib.pyplot as plt

def extract_and_plot_all_loops(filepath, zero_tol=None, min_points=50):
    """
    Extract each complete hysteresis loop (origin -> peak -> origin) from
    bidirectional cyclic loading data, and plot them individually.

    Parameters:
        filepath:   path to data file (two columns: displacement, force)
        zero_tol:   zero-crossing tolerance, default = 2% of max |displacement|
        min_points: minimum number of data points per loop (filters noise)
    Returns:
        loops: list of dicts with keys 'id','x','y','peak_disp','direction','energy'
    """
    # ---- 1. Load data ----
    data = np.loadtxt(filepath)
    x, y = data[:, 0], data[:, 1]
    print(f"Total points: {len(x)}, Displacement range: [{x.min():.2f}, {x.max():.2f}]")

    if zero_tol is None:
        zero_tol = 0.02 * np.max(np.abs(x))

    # ---- 2. Find zero-crossing points ----
    near_zero = np.abs(x) <= zero_tol
    diff_nz = np.diff(near_zero.astype(int))
    enters = np.where(diff_nz == 1)[0] + 1
    leaves = np.where(diff_nz == -1)[0] + 1

    if near_zero[0]:
        enters = np.concatenate([[0], enters])
    if near_zero[-1]:
        leaves = np.concatenate([leaves, [len(x)]])

    n = min(len(enters), len(leaves))
    zero_pts = sorted(set((enters[:n] + leaves[:n]) // 2))

    if not zero_pts or zero_pts[0] > min_points:
        zero_pts.insert(0, 0)

    # ---- 3. Extract loops between adjacent zero-crossings ----
    loops = []
    for i in range(len(zero_pts) - 1):
        s, e = zero_pts[i], zero_pts[i + 1]
        if e - s < min_points:
            continue
        xs, ys = x[s:e + 1], y[s:e + 1]
        peak = xs[np.argmax(np.abs(xs))]
        loops.append({
            'id': len(loops) + 1,
            'x': xs.copy(), 'y': ys.copy(),
            'peak_disp': peak,
            'direction': 'positive' if peak > 0 else 'negative',
            'energy': abs(np.trapz(ys, xs))
        })

    print(f"Extracted {len(loops)} complete hysteresis loops\n")

    # ---- 4. Overview plot ----
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, y, 'gray', linewidth=0.3, alpha=0.3, label='Full data')
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(loops), 1)))
    for i, lp in enumerate(loops):
        ax.plot(lp['x'], lp['y'], color=colors[i % 20], linewidth=1.2,
                alpha=0.8, label=f"#{lp['id']}")
    ax.set_xlabel('Displacement', fontsize=12)
    ax.set_ylabel('Force', fontsize=12)
    ax.set_title(f'All Hysteresis Loops Overview ({len(loops)} loops)', fontsize=14)
    if len(loops) <= 20:
        ax.legend(fontsize=7, ncol=3, loc='best')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    # ---- 5. Plot each loop individually ----
    n_loops = len(loops)
    n_cols = min(4, n_loops)
    n_rows = int(np.ceil(n_loops / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 4 * n_rows))
    axes = np.atleast_2d(axes)

    for i, lp in enumerate(loops):
        row, col = divmod(i, n_cols)
        ax = axes[row, col]
        c = 'red' if lp['direction'] == 'positive' else 'blue'

        ax.plot(lp['x'], lp['y'], color=c, linewidth=1.5)
        ax.fill(lp['x'], lp['y'], color=c, alpha=0.08)
        ax.plot(lp['x'][0], lp['y'][0], 'go', ms=6, zorder=5)

        pk = np.argmax(np.abs(lp['x']))
        ax.plot(lp['x'][pk], lp['y'][pk], 'k^', ms=8, zorder=5)

        ax.set_title(f"#{lp['id']}  {lp['direction']}\n"
                     f"peak={lp['peak_disp']:.1f}  E={lp['energy']:.2f}",
                     fontsize=9)
        ax.set_xlabel('Displacement', fontsize=8)
        ax.set_ylabel('Force', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
        ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')

    for i in range(n_loops, n_rows * n_cols):
        row, col = divmod(i, n_cols)
        axes[row, col].set_visible(False)

    plt.suptitle('Individual Hysteresis Loops', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.show()

    # ---- 6. Print summary ----
    print(f"{'No.':<5} {'Direction':<12} {'Peak Disp.':<14} {'Energy':<14} {'Points':<8}")
    print("-" * 55)
    for lp in loops:
        print(f"{lp['id']:<5} {lp['direction']:<12} {lp['peak_disp']:<14.4f} "
              f"{lp['energy']:<14.4f} {len(lp['x']):<8}")

    return loops

# ======================== Run ========================
if __name__ == "__main__":
    loops = extract_and_plot_all_loops("cableBearingData.txt", zero_tol=None, min_points=50)