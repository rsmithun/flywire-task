from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


DATASETS = {
    "BANC": "data/banc_626_edge_list.csv",
    "FAFB": "data/fafb_783_edge_list.csv",
    "MANC": "data/manc_1.2.1_edge_list.csv",
    "MAOL": "data/maol_1.1_edge_list.csv",
    "MCNS": "data/mcns_0.9_edge_list.csv",
}

TARGET_DATASETS = ["MANC", "MAOL", "MCNS"]


def load_nodes(path: Path) -> set[str]:
    nodes: set[str] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            source = row["source neuron id"]
            target = row["target neuron id"]
            nodes.add(source)
            nodes.add(target)
    return nodes


def load_edges_within(path: Path, allowed_nodes: set[str]) -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            source = row["source neuron id"]
            target = row["target neuron id"]
            if source in allowed_nodes and target in allowed_nodes:
                edges.add((source, target))
    return edges


def greedy_connected_size(root: int, common_adj: list[set[int]], conflict_adj: list[set[int]]) -> int:
    selected = {root}
    frontier = set(common_adj[root])

    while True:
        candidates = [
            vertex
            for vertex in frontier
            if vertex not in selected and not (conflict_adj[vertex] & selected)
        ]
        if not candidates:
            break

        candidates.sort(
            key=lambda vertex: (
                len(common_adj[vertex] - selected),
                len(common_adj[vertex]),
                -len(conflict_adj[vertex] & selected),
            ),
            reverse=True,
        )
        chosen = candidates[0]
        selected.add(chosen)
        frontier |= common_adj[chosen]
        frontier -= selected

    return len(selected)


def solve_rooted_milp(
    root: int,
    comp_nodes: list[str],
    common_edges: list[tuple[int, int]],
    conflict_pairs: list[tuple[int, int]],
) -> list[str]:
    n = len(comp_nodes)
    m = len(common_edges)
    max_flow = n - 1

    ix_x = 0
    ix_f = n
    num_vars = n + (2 * m)

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    lower: list[float] = []
    upper: list[float] = []

    # Fix the root to be selected.
    rows.append(len(lower))
    cols.append(ix_x + root)
    vals.append(1.0)
    lower.append(1.0)
    upper.append(1.0)

    # Pairwise conflicts: at most one endpoint may be selected.
    for a, b in conflict_pairs:
        row = len(lower)
        rows.extend([row, row])
        cols.extend([ix_x + a, ix_x + b])
        vals.extend([1.0, 1.0])
        lower.append(-np.inf)
        upper.append(1.0)

    incoming: list[list[int]] = [[] for _ in range(n)]
    outgoing: list[list[int]] = [[] for _ in range(n)]
    for edge_index, (u, v) in enumerate(common_edges):
        forward = ix_f + (2 * edge_index)
        reverse = ix_f + (2 * edge_index) + 1
        outgoing[u].append(forward)
        incoming[v].append(forward)
        outgoing[v].append(reverse)
        incoming[u].append(reverse)

    # One unit of flow reaches every selected non-root node.
    for vertex in range(n):
        row = len(lower)
        if vertex == root:
            for flow_index in outgoing[vertex]:
                rows.append(row)
                cols.append(flow_index)
                vals.append(1.0)
            for flow_index in incoming[vertex]:
                rows.append(row)
                cols.append(flow_index)
                vals.append(-1.0)
            for other in range(n):
                if other != root:
                    rows.append(row)
                    cols.append(ix_x + other)
                    vals.append(-1.0)
            lower.append(0.0)
            upper.append(0.0)
        else:
            for flow_index in incoming[vertex]:
                rows.append(row)
                cols.append(flow_index)
                vals.append(1.0)
            for flow_index in outgoing[vertex]:
                rows.append(row)
                cols.append(flow_index)
                vals.append(-1.0)
            rows.append(row)
            cols.append(ix_x + vertex)
            vals.append(-1.0)
            lower.append(0.0)
            upper.append(0.0)

    # Capacity constraints keep flow inside selected endpoints.
    for edge_index, (u, v) in enumerate(common_edges):
        forward = ix_f + (2 * edge_index)
        reverse = ix_f + (2 * edge_index) + 1

        for flow_index, a, b in ((forward, u, v), (reverse, v, u)):
            row = len(lower)
            rows.extend([row, row])
            cols.extend([flow_index, ix_x + a])
            vals.extend([1.0, -max_flow])
            lower.append(-np.inf)
            upper.append(0.0)

            row = len(lower)
            rows.extend([row, row])
            cols.extend([flow_index, ix_x + b])
            vals.extend([1.0, -max_flow])
            lower.append(-np.inf)
            upper.append(0.0)

    matrix = coo_matrix((vals, (rows, cols)), shape=(len(lower), num_vars)).tocsr()
    objective = np.zeros(num_vars)
    objective[:n] = -1.0

    integrality = np.zeros(num_vars, dtype=int)
    integrality[:n] = 1

    lower_bounds = np.zeros(num_vars)
    upper_bounds = np.full(num_vars, np.inf)
    upper_bounds[:n] = 1.0
    upper_bounds[ix_f:] = max_flow

    result = milp(
        c=objective,
        integrality=integrality,
        bounds=Bounds(lower_bounds, upper_bounds),
        constraints=LinearConstraint(matrix, np.array(lower), np.array(upper)),
        options={"disp": False, "time_limit": 120},
    )
    if result.x is None:
        raise RuntimeError(f"MILP failed for root {comp_nodes[root]}")

    selected_mask = np.rint(result.x[:n]).astype(int)
    return [comp_nodes[index] for index, value in enumerate(selected_mask) if value]


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce the shared 3-dataset circuit.")
    parser.add_argument("--data-root", type=Path, default=Path("."), help="Repository root containing data/")
    parser.add_argument("--out", type=Path, default=Path("network.csv"), help="Output CSV path")
    args = parser.parse_args()

    dataset_paths = {name: args.data_root / DATASETS[name] for name in TARGET_DATASETS}

    node_sets: dict[str, set[str]] = {}
    edge_sets: dict[str, set[tuple[str, str]]] = {}
    for name, path in dataset_paths.items():
        if not path.exists():
            raise FileNotFoundError(path)
        node_sets[name] = load_nodes(path)

    common_nodes = sorted(
        set.intersection(*(node_sets[name] for name in TARGET_DATASETS)),
        key=lambda value: int(value),
    )
    node_index = {node_id: index for index, node_id in enumerate(common_nodes)}

    allowed_nodes = set(common_nodes)
    for name, path in dataset_paths.items():
        edge_sets[name] = load_edges_within(path, allowed_nodes)

    common_directed = set.intersection(*(edge_sets[name] for name in TARGET_DATASETS))
    common_directed = {
        (source, target)
        for source, target in common_directed
        if source in node_index and target in node_index
    }
    undirected_common = {
        tuple(sorted((source, target), key=lambda value: int(value)))
        for source, target in common_directed
        if source != target
    }

    common_adj = [set() for _ in common_nodes]
    for source, target in undirected_common:
        u = node_index[source]
        v = node_index[target]
        common_adj[u].add(v)
        common_adj[v].add(u)

    seen = [False] * len(common_nodes)
    connected_components: list[list[int]] = []
    for start in range(len(common_nodes)):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        component: list[int] = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in common_adj[vertex]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(neighbor)
        connected_components.append(component)
    connected_components.sort(key=len, reverse=True)
    component = connected_components[0]
    component_set = set(component)
    comp_nodes = [common_nodes[index] for index in component]
    comp_index = {node_id: index for index, node_id in enumerate(comp_nodes)}

    component_adj = [set() for _ in comp_nodes]
    for source, target in undirected_common:
        if node_index[source] in component_set and node_index[target] in component_set:
            u = comp_index[source]
            v = comp_index[target]
            component_adj[u].add(v)
            component_adj[v].add(u)

    component_edges = [
        (comp_index[source], comp_index[target])
        for source, target in undirected_common
        if node_index[source] in component_set and node_index[target] in component_set
    ]

    # Build conflict pairs from any ordered pair whose presence pattern is not identical
    # across the three selected datasets.
    conflict_masks: defaultdict[tuple[int, int], int] = defaultdict(int)
    for bit, name in enumerate(TARGET_DATASETS):
        for source, target in edge_sets[name]:
            if source in node_index and target in node_index:
                u = node_index[source]
                v = node_index[target]
                if u in component_set and v in component_set:
                    conflict_masks[(u, v)] |= 1 << bit

    conflict_pairs = {
        tuple(sorted((comp_index[common_nodes[u]], comp_index[common_nodes[v]])))
        for (u, v), mask in conflict_masks.items()
        if mask not in (0, 7) and u != v
    }
    conflict_pairs_list = sorted(conflict_pairs)

    conflict_adj = [set() for _ in comp_nodes]
    for a, b in conflict_pairs_list:
        conflict_adj[a].add(b)
        conflict_adj[b].add(a)

    greedy_sizes = [
        (greedy_connected_size(root, component_adj, conflict_adj), root)
        for root in range(len(comp_nodes))
    ]
    greedy_sizes.sort(reverse=True)
    best_root = greedy_sizes[0][1]

    selected = solve_rooted_milp(best_root, comp_nodes, component_edges, conflict_pairs_list)
    if len(selected) != 12:
        raise RuntimeError(f"Expected 12 selected nodes, found {len(selected)}")

    selected_set = set(selected)
    induced_edges = []
    for name in TARGET_DATASETS:
        edges = {
            (source, target)
            for source, target in edge_sets[name]
            if source in selected_set and target in selected_set
        }
        induced_edges.append(edges)

    if not (induced_edges[0] == induced_edges[1] == induced_edges[2]):
        raise RuntimeError("Selected neurons do not induce identical directed subgraphs across the three datasets.")

    graph = nx.Graph()
    graph.add_nodes_from(selected)
    graph.add_edges_from({tuple(sorted(edge, key=lambda value: int(value))) for edge in induced_edges[0]})
    if not nx.is_connected(graph):
        raise RuntimeError("Selected neurons are not weakly connected.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(TARGET_DATASETS)
        for node_id in selected:
            writer.writerow([node_id, node_id, node_id])

    print(f"selected {len(selected)} neurons")
    print(f"root used for exact solve: {comp_nodes[best_root]}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()