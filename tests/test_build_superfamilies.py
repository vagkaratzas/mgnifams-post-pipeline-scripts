import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "build_superfamilies.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_superfamilies", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_filters_multiqc_csv_and_drops_row_column(tmp_path):
    module = load_module()
    input_csv = tmp_path / "similarity_mqc.csv"
    output_csv = tmp_path / "filtered.csv"
    input_csv.write_text(
        '# id: "family_similarities"\n'
        '# format: "csv"\n'
        "Row,Family Id 1,Family Id 2,Jaccard Score\n"
        "1,1,2,0.7\n"
        "2,1_loc,3,0.8\n"
        "3,4,5_loc,0.9\n"
        "4,2,3,0.4\n"
    )

    edges = module.read_filtered_edges(input_csv)
    module.write_filtered_edges(edges, output_csv)

    assert edges == [
        module.Edge("1", "2", 0.7),
        module.Edge("2", "3", 0.4),
    ]
    assert output_csv.read_text() == (
        "Family Id 1,Family Id 2,Jaccard Score\n"
        "1,2,0.7\n"
        "2,3,0.4\n"
    )


def test_builds_components_and_selects_weighted_representatives():
    module = load_module()
    edges = [
        module.Edge("1", "2", 0.6),
        module.Edge("2", "3", 0.9),
        module.Edge("10", "11", 0.5),
    ]

    graph = module.build_graph(edges)
    clusters = module.summarise_clusters(graph)

    assert clusters == [
        module.ClusterSummary("SF_1", "2", 3, ["1", "2", "3"]),
        module.ClusterSummary("SF_2", "10", 2, ["10", "11"]),
    ]


def test_appends_missing_ids_as_singletons_and_sorts_by_size(tmp_path):
    module = load_module()
    graph = module.build_graph(
        [
            module.Edge("1", "2", 0.6),
            module.Edge("2", "3", 0.9),
            module.Edge("10", "11", 0.5),
        ]
    )
    report_csv = tmp_path / "superfamily_statistics.csv"
    singleton_ids = tmp_path / "singleton_ids.txt"

    graph_clusters = module.summarise_clusters(graph)
    singletons = module.singleton_clusters_for_missing_ids(graph, min_id=1, max_id=12)
    clusters = module.sort_clusters_by_size(graph_clusters + singletons)
    module.write_singleton_ids(singletons, singleton_ids)
    module.write_cluster_report(clusters, report_csv)

    assert [cluster.family_rep_id for cluster in singletons] == [
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "12",
    ]
    assert singleton_ids.read_text() == "4\n5\n6\n7\n8\n9\n12\n"
    assert report_csv.read_text().splitlines()[:4] == [
        "Cluster Id,Family Rep Id,Family Size,Family Ids",
        "SF_1,2,3,1;2;3",
        "SF_2,10,2,10;11",
        "Singleton_4,4,1,4",
    ]


def test_groups_superfamily_sizes_into_ten_wide_bins():
    module = load_module()
    clusters = [
        module.ClusterSummary("SF_1", "1", 1, ["1"]),
        module.ClusterSummary("SF_2", "2", 2, ["2"]),
        module.ClusterSummary("SF_2", "2", 10, ["2"]),
        module.ClusterSummary("SF_3", "3", 11, ["3"]),
        module.ClusterSummary("SF_4", "4", 20, ["4"]),
        module.ClusterSummary("SF_5", "5", 21, ["5"]),
    ]

    assert module.low_size_distribution_bins(clusters) == [
        ("2", 1),
        ("10", 1),
    ]
    assert module.high_size_distribution_bins(clusters) == [
        ("11-20", 2),
        ("21-30", 1),
    ]


def test_writes_total_family_count(tmp_path):
    module = load_module()
    output_txt = tmp_path / "total_families.txt"
    clusters = [
        module.ClusterSummary("SF_1", "1", 3, ["1", "2", "3"]),
        module.ClusterSummary("Singleton_4", "4", 1, ["4"]),
        module.ClusterSummary("Singleton_5", "5", 1, ["5"]),
    ]

    total = module.write_total_family_count(clusters, output_txt)

    assert total == 5
    assert output_txt.read_text() == "5\n"


def test_default_output_path_uses_network_subdirectory(tmp_path):
    module = load_module()
    input_csv = tmp_path / "similarity_mqc.csv"

    assert module.default_output_path(input_csv, "superfamily_statistics.csv") == (
        tmp_path / "network" / "superfamily_statistics.csv"
    )
