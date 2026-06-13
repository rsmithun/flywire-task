from __future__ import annotations

import argparse
import csv
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.art3d import Line3DCollection


TARGET_DATASETS = ["MANC", "MAOL", "MCNS"]
EDGE_FILES = {
    "MANC": "data/manc_1.2.1_edge_list.csv",
    "MAOL": "data/maol_1.1_edge_list.csv",
    "MCNS": "data/mcns_0.9_edge_list.csv",
}
MANC_META_URL = (
    "https://storage.googleapis.com/"
    "lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/manc_121/manc_121_meta.feather"
)
MANC_SWC_BASE = (
    "https://storage.googleapis.com/"
    "lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/manc_121/manc_manc_space_swc"
)


def download_file(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        urllib.request.urlretrieve(url, path)


def load_selected_ids(network_csv: Path) -> list[str]:
    with network_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [row["MANC"] for row in reader]


def load_annotations(cache_dir: Path, selected: list[str]) -> dict[str, dict[str, str]]:
    feather_path = cache_dir / "manc_121_meta.feather"
    download_file(MANC_META_URL, feather_path)
    annotations = pd.read_feather(feather_path)
    annotations["manc_121_id"] = annotations["manc_121_id"].astype(str)
    subset = annotations[annotations["manc_121_id"].isin(selected)].copy()

    rows: dict[str, dict[str, str]] = {}
    for _, row in subset.iterrows():
        rows[str(row["manc_121_id"])] = {
            "cell_type": str(row.get("cell_type", "")),
            "super_class": str(row.get("super_class", "")),
            "nt": str(row.get("neurotransmitter_predicted", "")).lower(),
        }
    return rows


def load_common_edges(selected: set[str], data_root: Path) -> set[tuple[str, str]]:
    common: set[tuple[str, str]] | None = None
    for name in TARGET_DATASETS:
        path = data_root / EDGE_FILES[name]
        edges: set[tuple[str, str]] = set()
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                source = row["source neuron id"]
                target = row["target neuron id"]
                if source in selected and target in selected:
                    edges.add((source, target))
        common = edges if common is None else common & edges
    return common or set()


def plot_network(selected: list[str], annotations: dict[str, dict[str, str]], edges: set[tuple[str, str]], out_path: Path) -> None:
    palette = {
        "acetylcholine": "#f28e2b",
        "gaba": "#4e79a7",
        "glutamate": "#59a14f",
        "sensory": "#e15759",
    }

    positions = {
        "10114": (0.0, 0.0),
        "37836": (1.6, 0.25),
        "28795": (-1.3, 0.95),
        "34613": (-2.2, -0.2),
        "11595": (-1.4, -1.15),
        "21042": (-2.2, -1.35),
        "24996": (0.95, 1.2),
        "33667": (1.55, 0.95),
        "22626": (2.05, 0.55),
        "19404": (2.2, -0.05),
        "21401": (2.0, -0.7),
        "15447": (0.85, -1.1),
    }

    fig, ax = plt.subplots(figsize=(11, 8.5), dpi=220)
    ax.set_facecolor("white")

    for source, target in edges:
        curvature = 0.0
        if source == "10114" and target == "37836":
            curvature = 0.22
        elif source == "37836" and target == "10114":
            curvature = -0.22

        arrow = FancyArrowPatch(
            posA=positions[source],
            posB=positions[target],
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.6,
            color="#6b7280",
            alpha=0.85,
            connectionstyle=f"arc3,rad={curvature}",
        )
        ax.add_patch(arrow)

    for node_id in selected:
        annotation = annotations[node_id]
        nt = annotation["nt"]
        color = palette.get(nt, palette["sensory"] if annotation["super_class"] == "sensory" else "#999999")
        size = 1500 if node_id == "10114" else 900 if node_id in {"34613", "37836", "28795"} else 760
        edgecolor = "#111827" if node_id == "10114" else "white"
        ax.scatter(*positions[node_id], s=size, c=color, edgecolors=edgecolor, linewidths=1.8, zorder=3)
        ax.text(
            positions[node_id][0],
            positions[node_id][1],
            f"{node_id}\n{annotation['cell_type']}",
            ha="center",
            va="center",
            fontsize=8.2,
            color="black",
            zorder=4,
        )

    ax.text(
        0.02,
        0.98,
        "Exact 12-node shared circuit across MANC / MAOL / MCNS",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.925,
        "Central hub: 10114 / vMS17 (GABAergic bilateral interconnecting interneuron)",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        color="#374151",
    )

    handles = [
        Line2D([0], [0], color=color, lw=2.5, label=label)
        for label, color in [
            ("Sensory afferent", palette["sensory"]),
            ("Acetylcholine", palette["acetylcholine"]),
            ("GABA", palette["gaba"]),
            ("Glutamate", palette["glutamate"]),
        ]
    ]
    ax.legend(handles=handles, loc="upper right", frameon=True, framealpha=0.95, fontsize=9, title="Node color")

    ax.set_xlim(-2.9, 2.9)
    ax.set_ylim(-1.9, 1.75)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def load_swc(path: Path) -> tuple[list[list[tuple[float, float, float]]], list[tuple[float, float, float]], tuple[float, float, float]]:
    points: dict[int, tuple[float, float, float]] = {}
    segments: list[list[tuple[float, float, float]]] = []
    all_points: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            point_id = int(parts[0])
            x, y, z = map(float, parts[2:5])
            parent = int(parts[6])
            current = (x, y, z)
            points[point_id] = current
            all_points.append(current)
            if parent != -1 and parent in points:
                segments.append([points[parent], current])
    root_id = min(points.keys())
    return segments, all_points, points[root_id]


def plot_skeletons(selected: list[str], annotations: dict[str, dict[str, str]], cache_dir: Path, out_path: Path) -> None:
    palette = {
        "acetylcholine": "#f28e2b",
        "gaba": "#4e79a7",
        "glutamate": "#59a14f",
        "sensory": "#e15759",
    }

    swc_dir = cache_dir / "manc_swc"
    swc_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(11, 9), dpi=220)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")

    all_points: list[tuple[float, float, float]] = []
    for node_id in selected:
        swc_path = swc_dir / f"{node_id}.swc"
        download_file(f"{MANC_SWC_BASE}/{node_id}.swc", swc_path)
        segments, points, root_point = load_swc(swc_path)
        annotation = annotations[node_id]
        nt = annotation["nt"]
        color = palette.get(nt, palette["sensory"] if annotation["super_class"] == "sensory" else "#999999")
        if segments:
            ax.add_collection3d(Line3DCollection(segments, colors=color, linewidths=0.25, alpha=0.12))
        ax.scatter([root_point[0]], [root_point[1]], [root_point[2]], s=18 if node_id != "10114" else 30, color=color, alpha=0.95)
        ax.text(root_point[0], root_point[1], root_point[2], node_id, fontsize=6.5, color="#111827")
        all_points.extend(points)

    coords = np.array(all_points)
    lower = coords.min(axis=0)
    upper = coords.max(axis=0)
    center = (lower + upper) / 2.0
    radius = (upper - lower).max() / 2.0
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.view_init(elev=20, azim=-58)
    ax.set_axis_off()
    ax.set_title("Native MANC SWC skeletons for the shared 12-node circuit", fontsize=15, pad=14)

    handles = [
        Line2D([0], [0], color=color, lw=2.5, label=label)
        for label, color in [
            ("Sensory afferent", palette["sensory"]),
            ("Acetylcholine", palette["acetylcholine"]),
            ("GABA", palette["gaba"]),
            ("Glutamate", palette["glutamate"]),
        ]
    ]
    ax.legend(handles=handles, loc="upper left", frameon=True, framealpha=0.95, fontsize=9)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerate the published circuit figures.")
    parser.add_argument("--data-root", type=Path, default=Path("."), help="Repository root containing data/")
    parser.add_argument("--network", type=Path, default=Path("network.csv"), help="Matched-neuron CSV")
    parser.add_argument("--out", type=Path, default=Path("figures"), help="Output directory")
    args = parser.parse_args()

    selected = load_selected_ids(args.network)
    selected_set = set(selected)
    cache_dir = args.out / ".cache"
    annotations = load_annotations(cache_dir, selected)
    edges = load_common_edges(selected_set, args.data_root)

    plot_network(selected, annotations, edges, args.out / "network.png")
    plot_skeletons(selected, annotations, cache_dir, args.out / "manc_skeletons_3d.png")

    print(f"wrote {args.out / 'network.png'}")
    print(f"wrote {args.out / 'manc_skeletons_3d.png'}")


if __name__ == "__main__":
    main()