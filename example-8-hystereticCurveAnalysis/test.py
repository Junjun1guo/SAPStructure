import numpy as np
import matplotlib.pyplot as plt

def extract_skeleton_curve(filepath, zero_tol=None, min_points=50):
    """
    Extract the skeleton (backbone) curve from hysteresis data and plot it
    over the original hysteresis loops.

    The skeleton curve connects the peak force-displacement points
    from each half-cycle of the hysteresis response.

    Parameters:
        filepath:   data file path (two columns: displacement, force)
        zero_tol:   zero-crossing tolerance, default = 2% of max |displacement|
        min_points: minimum data points per half-cycle to be considered valid
    Returns:
        skeleton_pos: ndarray (N, 2) - positive branch [disp, force]
        skeleton_neg: ndarray (M, 2) - negative branch [disp, force]
    """
    # ---- 1. Load data ----
    data = np.loadtxt(filepath)
    x, y = data[:, 0], data[:, 1]
    print(f"Total points: {len(x)}, Disp range: [{x.min():.2f}, {x.max():.2f}]")

    if zero_tol is None:
        zero_tol = 0.02 * np.max(np.abs(x))

    # ---- 2. Find zero-crossing regions ----
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

    # ---- 3. Extract peak points from each half-cycle ----
    pos_peaks_disp = []
    pos_peaks_force = []
    neg_peaks_disp = []
    neg_peaks_force = []

    for i in range(len(zero_pts) - 1):
        s, e = zero_pts[i], zero_pts[i + 1]
        if e - s < min_points:
            continue

        xs = x[s:e + 1]
        ys = y[s:e + 1]

        # Find the index of the peak absolute displacement in this half-cycle
        idx_peak = np.argmax(np.abs(xs))
        peak_disp = xs[idx_peak]
        peak_force = ys[idx_peak]

        if peak_disp > 0:
            pos_peaks_disp.append(peak_disp)
            pos_peaks_force.append(peak_force)
        else:
            neg_peaks_disp.append(peak_disp)
            neg_peaks_force.append(peak_force)

    # ---- 4. Build skeleton curve arrays ----
    # Add origin point (0, 0) to both branches for a complete skeleton
    # Positive branch: sort by displacement ascending
    pos_peaks_disp = np.array(pos_peaks_disp)
    pos_peaks_force = np.array(pos_peaks_force)
    neg_peaks_disp = np.array(neg_peaks_disp)
    neg_peaks_force = np.array(neg_peaks_force)

    # Sort positive branch by displacement (ascending)
    if len(pos_peaks_disp) > 0:
        sort_pos = np.argsort(pos_peaks_disp)
        pos_peaks_disp = pos_peaks_disp[sort_pos]
        pos_peaks_force = pos_peaks_force[sort_pos]

        # Remove duplicate displacement levels: keep the one with max force
        unique_pos_disp = []
        unique_pos_force = []
        visited = set()
        # Round to a tolerance to group similar displacement amplitudes
        disp_round = np.round(pos_peaks_disp, decimals=0)
        for d_val in np.unique(disp_round):
            mask = disp_round == d_val
            idx_max_f = np.argmax(pos_peaks_force[mask])
            unique_pos_disp.append(pos_peaks_disp[mask][idx_max_f])
            unique_pos_force.append(pos_peaks_force[mask][idx_max_f])
        pos_peaks_disp = np.array(unique_pos_disp)
        pos_peaks_force = np.array(unique_pos_force)

    # Sort negative branch by displacement (descending, i.e., from 0 toward most negative)
    if len(neg_peaks_disp) > 0:
        sort_neg = np.argsort(neg_peaks_disp)[::-1]
        neg_peaks_disp = neg_peaks_disp[sort_neg]
        neg_peaks_force = neg_peaks_force[sort_neg]

        # Remove duplicate displacement levels: keep the one with min (most negative) force
        unique_neg_disp = []
        unique_neg_force = []
        disp_round = np.round(neg_peaks_disp, decimals=0)
        for d_val in np.unique(disp_round)[::-1]:
            mask = disp_round == d_val
            idx_min_f = np.argmin(neg_peaks_force[mask])
            unique_neg_disp.append(neg_peaks_disp[mask][idx_min_f])
            unique_neg_force.append(neg_peaks_force[mask][idx_min_f])
        neg_peaks_disp = np.array(unique_neg_disp)
        neg_peaks_force = np.array(unique_neg_force)

    # Prepend origin to both branches
    skeleton_pos = np.column_stack([
        np.concatenate([[0.0], pos_peaks_disp]),
        np.concatenate([[0.0], pos_peaks_force])
    ])
    skeleton_neg = np.column_stack([
        np.concatenate([[0.0], neg_peaks_disp]),
        np.concatenate([[0.0], neg_peaks_force])
    ])

    print(f"Skeleton curve: {len(skeleton_pos)} positive points, {len(skeleton_neg)} negative points")

    # ---- 5. Plot ----
    fig, ax = plt.subplots(figsize=(12, 7))

    # Raw hysteresis in gray
    ax.plot(x, y, color='gray', linewidth=0.3, alpha=0.5, label='Hysteresis Curve')

    # Skeleton curve in blue
    ax.plot(skeleton_pos[:, 0], skeleton_pos[:, 1],
            'b-o', linewidth=2.0, markersize=6, markerfacecolor='blue',
            markeredgecolor='white', markeredgewidth=0.8, zorder=5,
            label='Skeleton Curve (Positive)')
    ax.plot(skeleton_neg[:, 0], skeleton_neg[:, 1],
            'b-s', linewidth=2.0, markersize=6, markerfacecolor='blue',
            markeredgecolor='white', markeredgewidth=0.8, zorder=5,
            label='Skeleton Curve (Negative)')

    # Mark peak points
    ax.plot(skeleton_pos[1:, 0], skeleton_pos[1:, 1], 'ro', ms=4, zorder=6, alpha=0.6)
    ax.plot(skeleton_neg[1:, 0], skeleton_neg[1:, 1], 'ro', ms=4, zorder=6, alpha=0.6)

    # Annotate peak values
    for i in range(1, len(skeleton_pos)):
        ax.annotate(f'({skeleton_pos[i, 0]:.1f}, {skeleton_pos[i, 1]:.1f})',
                    xy=(skeleton_pos[i, 0], skeleton_pos[i, 1]),
                    xytext=(5, 8), textcoords='offset points',
                    fontsize=7, color='blue', alpha=0.8)
    for i in range(1, len(skeleton_neg)):
        ax.annotate(f'({skeleton_neg[i, 0]:.1f}, {skeleton_neg[i, 1]:.1f})',
                    xy=(skeleton_neg[i, 0], skeleton_neg[i, 1]),
                    xytext=(5, -12), textcoords='offset points',
                    fontsize=7, color='blue', alpha=0.8)

    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
    ax.set_xlabel('Displacement', fontsize=12)
    ax.set_ylabel('Force', fontsize=12)
    ax.set_title('Hysteresis Curve with Skeleton (Backbone) Curve', fontsize=14)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    # ---- 6. Print skeleton data ----
    print("\n--- Positive Branch ---")
    print(f"{'Displacement':<15} {'Force':<15}")
    print("-" * 30)
    for i in range(len(skeleton_pos)):
        print(f"{skeleton_pos[i, 0]:<15.4f} {skeleton_pos[i, 1]:<15.4f}")

    print("\n--- Negative Branch ---")
    print(f"{'Displacement':<15} {'Force':<15}")
    print("-" * 30)
    for i in range(len(skeleton_neg)):
        print(f"{skeleton_neg[i, 0]:<15.4f} {skeleton_neg[i, 1]:<15.4f}")

    return skeleton_pos, skeleton_neg

# ======================== Usage ========================
if __name__ == "__main__":
    skeleton_pos, skeleton_neg = extract_skeleton_curve("cableBearingData.txt")