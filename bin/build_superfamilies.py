#!/usr/bin/env python3

import argparse
import csv
import os
from collections import defaultdict, deque
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Set, Tuple


DEFAULT_INPUT = (
    "assets/mgnifams_v2_results/generate_families/similarity_mqc.csv"
)
DEFAULT_MIN_FAMILY_ID = 1
DEFAULT_MAX_FAMILY_ID = 35459


class Edge(NamedTuple):
    family_id_1: str
    family_id_2: str
    weight: float


class ClusterSummary(NamedTuple):
    cluster_id: str
    family_rep_id: str
    family_size: int
    family_ids: List[str]


Graph = Dict[str, Dict[str, float]]


def family_sort_key(family_id: str) -> Tuple[int, object]:
    """Sort numeric MGnifam ids numerically, while still accepting text ids."""
    try:
        return (0, int(family_id))
    except ValueError:
        return (1, family_id)


def csv_data_lines(path):
    """Yield data rows from a MultiQC CSV, skipping leading metadata comments."""
    with open(path, newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            yield line


def read_filtered_edges(path) -> List[Edge]:
    """Read valid similarity edges, excluding localized family ids.

    The input is a MultiQC-flavoured CSV whose first column is a display-only
    row number. Edges where either endpoint contains "_" are discarded because
    those ids are not canonical family ids for this superfamily graph.
    """
    reader = csv.DictReader(csv_data_lines(path))
    required_columns = {"Family Id 1", "Family Id 2", "Jaccard Score"}
    missing_columns = required_columns.difference(reader.fieldnames or [])
    if missing_columns:
        raise ValueError(
            f"Missing required columns in {path}: {', '.join(sorted(missing_columns))}"
        )

    edges = []
    for row in reader:
        family_id_1 = row["Family Id 1"].strip()
        family_id_2 = row["Family Id 2"].strip()
        if "_" in family_id_1 or "_" in family_id_2:
            continue
        edges.append(Edge(family_id_1, family_id_2, float(row["Jaccard Score"])))
    return edges


def write_filtered_edges(edges: Iterable[Edge], output_csv) -> None:
    """Write the filtered edge list without the original MultiQC Row column."""
    with open(output_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Family Id 1", "Family Id 2", "Jaccard Score"])
        for family_id_1, family_id_2, weight in edges:
            writer.writerow([family_id_1, family_id_2, weight])


def build_graph(edges: Iterable[Edge]) -> Graph:
    """Build an undirected weighted graph from family similarity edges."""
    graph = defaultdict(dict)
    for family_id_1, family_id_2, weight in edges:
        graph[family_id_1][family_id_2] = weight
        graph[family_id_2][family_id_1] = weight
    return graph


def connected_components(graph: Graph) -> List[List[str]]:
    """Return connected components sorted by their smallest family id."""
    visited: Set[str] = set()
    components = []

    for start in sorted(graph, key=family_sort_key):
        if start in visited:
            continue
        component = []
        queue = deque([start])
        visited.add(start)

        while queue:
            family_id = queue.popleft()
            component.append(family_id)
            for neighbour in sorted(graph[family_id], key=family_sort_key):
                if neighbour not in visited:
                    visited.add(neighbour)
                    queue.append(neighbour)

        components.append(sorted(component, key=family_sort_key))

    return components


def representative_for_component(graph: Graph, family_ids: List[str]) -> str:
    """Select the family nearest the component centroid by weighted degree.

    Without family-level embeddings or coordinates, the centroid is represented
    by the most strongly connected node: the family with the largest sum of
    Jaccard weights to other families in the same connected component.
    """
    family_set = set(family_ids)

    def weighted_degree(family_id):
        return sum(
            weight
            for neighbour, weight in graph[family_id].items()
            if neighbour in family_set
        )

    return sorted(
        family_ids,
        key=lambda family_id: (-weighted_degree(family_id), family_sort_key(family_id)),
    )[0]


def summarise_clusters(graph: Graph) -> List[ClusterSummary]:
    """Create superfamily summaries from the graph's connected components."""
    summaries = []
    for index, family_ids in enumerate(connected_components(graph), start=1):
        summaries.append(
            ClusterSummary(
                f"SF_{index}",
                representative_for_component(graph, family_ids),
                len(family_ids),
                family_ids,
            )
        )
    return summaries


def singleton_clusters_for_missing_ids(
    graph: Graph,
    min_id: int = DEFAULT_MIN_FAMILY_ID,
    max_id: int = DEFAULT_MAX_FAMILY_ID,
) -> List[ClusterSummary]:
    """Create singleton superfamilies for canonical ids absent from all edges.

    The MGnifam v2 canonical id space is fixed here as 1..35459. Any id in that
    range missing from both edge columns has no observed similarity edge, so it
    must still be represented in the aggregate superfamily report as size 1.
    """
    present_ids = set(graph)
    return [
        ClusterSummary(f"Singleton_{family_id}", str(family_id), 1, [str(family_id)])
        for family_id in range(min_id, max_id + 1)
        if str(family_id) not in present_ids
    ]


def sort_clusters_by_size(clusters: Iterable[ClusterSummary]) -> List[ClusterSummary]:
    """Sort report rows by descending superfamily size, then representative id."""
    return sorted(
        clusters,
        key=lambda cluster: (
            -cluster.family_size,
            family_sort_key(cluster.family_rep_id),
            cluster.cluster_id,
        ),
    )


def write_singleton_ids(singletons: Iterable[ClusterSummary], output_txt) -> None:
    """Write one missing canonical family id per line."""
    with open(output_txt, "w") as handle:
        for singleton in singletons:
            handle.write(f"{singleton.family_rep_id}\n")


def write_cluster_report(clusters: Iterable[ClusterSummary], output_csv) -> None:
    """Write per-superfamily representative and size statistics."""
    with open(output_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Cluster Id", "Family Rep Id", "Family Size", "Family Ids"])
        for cluster_id, family_rep_id, family_size, family_ids in clusters:
            writer.writerow(
                [cluster_id, family_rep_id, family_size, ";".join(family_ids)]
            )


def low_size_distribution_bins(clusters: Iterable[ClusterSummary]) -> List[Tuple[str, int]]:
    """Count non-singleton superfamilies with exact sizes from 2 through 10."""
    counts = defaultdict(int)
    for cluster in clusters:
        if 2 <= cluster.family_size <= 10:
            counts[cluster.family_size] += 1

    return [(str(size), counts[size]) for size in range(2, 11) if counts[size]]


def high_size_distribution_bins(clusters: Iterable[ClusterSummary]) -> List[Tuple[str, int]]:
    """Count superfamilies of size 11+ in 10-wide family-size bins."""
    counts = defaultdict(int)
    for cluster in clusters:
        if cluster.family_size < 11:
            continue
        lower = ((cluster.family_size - 1) // 10) * 10 + 1
        upper = lower + 9
        counts[(lower, upper)] += 1

    return [
        (f"{lower}-{upper}", counts[(lower, upper)])
        for lower, upper in sorted(counts)
    ]


def write_distribution_plot(
    binned_counts: List[Tuple[str, int]],
    output_png,
    xlabel: str,
    title: str,
) -> None:
    """Render a barplot from precomputed superfamily size bins."""
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [label for label, count in binned_counts]
    counts = [count for label, count in binned_counts]

    width = max(8, min(24, len(labels) * 0.35))
    fig, ax = plt.subplots(figsize=(width, 5))
    ax.bar(labels, counts, color="#3572a5")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of superfamilies")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=90 if len(labels) > 12 else 0)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200)
    plt.close(fig)


def write_size_distribution_plots(
    small_plot_clusters: List[ClusterSummary],
    large_plot_clusters: List[ClusterSummary],
    low_output_png,
    high_output_png,
) -> None:
    """Render exact 2-10 sizes separately from 11+ ranges.

    The first plot intentionally receives only graph-derived clusters so the
    large singleton population does not dominate the visible 2-10 distribution.
    """
    write_distribution_plot(
        low_size_distribution_bins(small_plot_clusters),
        low_output_png,
        "Superfamily size",
        "Superfamily size distribution: 2-10",
    )
    write_distribution_plot(
        high_size_distribution_bins(large_plot_clusters),
        high_output_png,
        "Superfamily size range",
        "Superfamily size distribution: 11+",
    )


def write_total_family_count(clusters: Iterable[ClusterSummary], output_txt) -> int:
    """Write the sum of Family Size from the final superfamily report."""
    total = sum(cluster.family_size for cluster in clusters)
    with open(output_txt, "w") as handle:
        handle.write(f"{total}\n")
    return total


def default_output_path(input_csv, filename):
    return input_csv.parent / "network" / filename


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build connected-component superfamilies from an MGnifam family "
            "similarity MultiQC CSV."
        )
    )
    parser.add_argument(
        "input_csv",
        nargs="?",
        type=Path,
        default=Path(DEFAULT_INPUT),
        help=f"Input similarity_mqc.csv file. Defaults to {DEFAULT_INPUT}.",
    )
    parser.add_argument(
        "--filtered-csv",
        help="Filtered edge-list CSV to write. Defaults next to input CSV.",
    )
    parser.add_argument(
        "--report-csv",
        help="Superfamily statistics CSV to write. Defaults next to input CSV.",
    )
    parser.add_argument(
        "--singleton-ids",
        help="Missing canonical family ids TXT to write. Defaults next to input CSV.",
    )
    parser.add_argument(
        "--plot-small-png",
        help="Exact 1-10 superfamily size distribution PNG. Defaults next to input CSV.",
    )
    parser.add_argument(
        "--plot-large-png",
        help="Range-binned 11+ superfamily size distribution PNG. Defaults next to input CSV.",
    )
    parser.add_argument(
        "--total-families",
        help="Total family count sanity-check TXT to write. Defaults next to input CSV.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_csv = args.input_csv
    filtered_csv = Path(
        args.filtered_csv
        or default_output_path(input_csv, "similarity_filtered_superfamilies.csv")
    )
    report_csv = Path(
        args.report_csv
        or default_output_path(input_csv, "superfamily_statistics.csv")
    )
    singleton_ids = Path(
        args.singleton_ids
        or default_output_path(input_csv, "singleton_ids.txt")
    )
    plot_small_png = Path(
        args.plot_small_png
        or default_output_path(input_csv, "superfamily_size_distribution_2_to_10.png")
    )
    plot_large_png = Path(
        args.plot_large_png
        or default_output_path(input_csv, "superfamily_size_distribution_11_plus.png")
    )
    total_families = Path(
        args.total_families
        or default_output_path(input_csv, "total_families.txt")
    )

    for output_path in (
        filtered_csv,
        report_csv,
        singleton_ids,
        plot_small_png,
        plot_large_png,
        total_families,
    ):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    edges = read_filtered_edges(input_csv)
    write_filtered_edges(edges, filtered_csv)

    graph = build_graph(edges)
    graph_clusters = summarise_clusters(graph)
    singleton_clusters = singleton_clusters_for_missing_ids(graph)
    clusters = sort_clusters_by_size(graph_clusters + singleton_clusters)
    write_singleton_ids(singleton_clusters, singleton_ids)
    write_cluster_report(clusters, report_csv)
    total_family_count = write_total_family_count(clusters, total_families)
    write_size_distribution_plots(
        graph_clusters,
        clusters,
        plot_small_png,
        plot_large_png,
    )

    print(f"Filtered edges: {len(edges)} -> {filtered_csv}")
    print(f"Singleton family ids: {len(singleton_clusters)} -> {singleton_ids}")
    print(f"Superfamilies: {len(clusters)} -> {report_csv}")
    print(f"Total families: {total_family_count} -> {total_families}")
    print(f"Size distribution plot 2-10: {plot_small_png}")
    print(f"Size distribution plot 11+: {plot_large_png}")


if __name__ == "__main__":
    main()
