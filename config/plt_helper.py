from matplotlib import rcParams
rcParams.update({'font.size': 20})

import os
import numpy as np
from matplotlib.transforms import Bbox

def save_each_axes(fig, axs, out_dir, prefix="subplot", dpi=300, pad_frac=0.02,
                   include_outside_legend=True):
    """
    Save each Axes from an existing figure by cropping the figure canvas
    to that Axes' (and legend's) tight bounding box.

    Works with lines, fill_between (PolyCollection), images, etc.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    axs : array-like of Axes
    out_dir : str
    prefix : str
    dpi : int
    pad_frac : float
        Extra fractional padding around the axes bbox.
    include_outside_legend : bool
        If True, expand bbox to include legend even if it's outside the axes.
    """
    os.makedirs(out_dir, exist_ok=True)
    axs = np.ravel(axs)

    # Ensure tightbbox values are up-to-date
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    for i, ax in enumerate(axs):
        bbox = ax.get_tightbbox(renderer)

        # Optionally include legend box if placed outside the axes
        leg = ax.get_legend()
        if include_outside_legend and leg is not None:
            leg_bbox = leg.get_window_extent(renderer)
            bbox = Bbox.union([bbox, leg_bbox])

        # Convert from pixels to inches and add a little padding
        bbox_in = bbox.transformed(fig.dpi_scale_trans.inverted())
        bbox_in = bbox_in.expanded(1 + pad_frac, 1 + pad_frac)

        out = os.path.join(out_dir, f"{prefix}_{i}.png")
        out_pdf = os.path.join(out_dir, f"{prefix}_{i}.pdf")
        fig.savefig(out, dpi=dpi, bbox_inches=bbox_in)
        fig.savefig(out_pdf, dpi=dpi, bbox_inches=bbox_in)
        print(f"Saved {out}")

def save_subplots(axes, plot_dir, equalr=False):
    # Save individual subplots
    for i, ax_orig in enumerate(axes):
        # Create a new figure and axis
        fig_single, ax_single = plt.subplots(figsize=(8, 6))

        # Copy lines from the original axis
        for line in ax_orig.get_lines():
            ax_single.plot(*line.get_data(),
                           label=line.get_label(),
                           color=line.get_color(),
                           linestyle=line.get_linestyle(),
                           marker=line.get_marker(),
                           markersize=line.get_markersize())

        # --- Copy fill_between patches ---
        for coll in ax_orig.collections:
            if isinstance(coll, PolyCollection):
                # Extract vertices of the polygon(s)
                for path in coll.get_paths():
                    verts = path.vertices
                    x = verts[:, 0]
                    y = verts[:, 1]
                    # Because fill_between creates top and bottom polygons,
                    # you can just re-plot them
                    ax_single.fill_between(
                        x, y.min(), y.max(),
                        color=coll.get_facecolor()[0],
                        alpha=coll.get_alpha()
                    )

        # Copy axis formatting
        ax_single.set_xscale(ax_orig.get_xscale())
        ax_single.set_yscale(ax_orig.get_yscale())
        ax_single.set_xlabel(ax_orig.get_xlabel())
        ax_single.set_ylabel(ax_orig.get_ylabel())
        ax_single.set_title(ax_orig.get_title())
        if equalr: 
            ax_single.set_aspect('equal', adjustable='box')
        ax_single.grid(True)

        # Copy legend
        handles, labels = ax_orig.get_legend_handles_labels()
        if handles:
            # ax_single.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            # ax_single.legend(handles, labels, loc='lower left', frameon=True)
            # indp.legend(loc='lower left', frameon=True)
            # ax_single.legend(handles, labels)
            ax_single.legend(loc='best', frameon=True)

        # Save the individual figure
        # fig_single.tight_layout()
        fig_single.savefig(f"{plot_dir}_{i}.pdf", dpi=300, bbox_inches='tight')
        fig_single.savefig(f"{plot_dir}_{i}.png", dpi=300, bbox_inches='tight')
        plt.close(fig_single)