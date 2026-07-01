import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "biome_analysis.py"


def load_module():
    """Import the standalone script without requiring `bin` to be a package."""
    spec = importlib.util.spec_from_file_location("biome_analysis", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parse_biome_blob_accepts_sqlite_bytes():
    module = load_module()
    # SQLite may return BLOB columns as bytes; the public parser should accept
    # that shape even though the analysis hot path uses record dictionaries.
    blob = (
        b"ids,labels,parents,counts\n"
        b"root,root,,2\n"
        b"root:Environmental,Environmental,root,2\n"
        b"root:Environmental:Soil,Soil,root:Environmental,2\n"
    )

    df = module.parse_biome_blob(blob)

    assert list(df["labels"]) == ["root", "Environmental", "Soil"]
    assert list(df["counts"]) == [2, 2, 2]


def test_analyse_biomes_handles_no_parseable_families():
    module = load_module()

    result = module.analyse_biomes([])

    assert result["total_families"] == 0
    assert result["leaf_df"].empty
    assert result["exclusive"] == {}


def test_leaf_distribution_counts_families_not_sequences():
    module = load_module()
    rows = [
        (
            "MGF0001",
            (
                "ids,labels,parents,counts\n"
                "root,root,,100\n"
                "root:Environmental,Environmental,root,100\n"
                "root:Environmental:Soil,Soil,root:Environmental,100\n"
            ),
        ),
        (
            "MGF0002",
            (
                "ids,labels,parents,counts\n"
                "root,root,,5\n"
                "root:Environmental,Environmental,root,5\n"
                "root:Environmental:Soil,Soil,root:Environmental,5\n"
            ),
        ),
    ]

    result = module.analyse_biomes(rows)
    soil_count = (
        result["leaf_df"]
        .set_index("label")
        .loc["root:Environmental:Soil", "count"]
    )

    assert soil_count == 2


def test_exclusive_biomes_are_leaf_biomes_seen_in_one_family():
    module = load_module()
    rows = [
        (
            "MGF0001",
            (
                "ids,labels,parents,counts\n"
                "root,root,,10\n"
                "root:Environmental,Environmental,root,10\n"
                "root:Environmental:Salt crystallizer pond,Salt crystallizer pond,root:Environmental,1\n"
                "root:Environmental:Soil,Soil,root:Environmental,9\n"
            ),
        ),
        (
            "MGF0002",
            (
                "ids,labels,parents,counts\n"
                "root,root,,5\n"
                "root:Environmental,Environmental,root,5\n"
                "root:Environmental:Soil,Soil,root:Environmental,5\n"
            ),
        ),
    ]

    result = module.analyse_biomes(rows)

    assert result["exclusive_leaf_biomes"] == {
        "root:Environmental:Salt crystallizer pond": ["MGF0001"]
    }


def test_leaf_distribution_keeps_duplicate_terminal_labels_as_separate_paths():
    module = load_module()
    # Regression guard: identical terminal labels from different branches must
    # not be collapsed into one Mixed biome.
    rows = [
        (
            "MGF0001",
            (
                "ids,labels,parents,counts\n"
                "root,root,,1\n"
                "root:Environmental,Environmental,root,1\n"
                "root:Environmental:Soil,Soil,root:Environmental,1\n"
            ),
        ),
        (
            "MGF0002",
            (
                "ids,labels,parents,counts\n"
                "root,root,,1\n"
                "root:Host-associated,Host-associated,root,1\n"
                "root:Host-associated:Plants,Plants,root:Host-associated,1\n"
                "root:Host-associated:Plants:Soil,Soil,root:Host-associated:Plants,1\n"
            ),
        ),
    ]

    result = module.analyse_biomes(rows)
    leaf_df = result["leaf_df"].set_index("label")

    assert set(leaf_df.index) == {
        "root:Environmental:Soil",
        "root:Host-associated:Plants:Soil",
    }
    assert leaf_df.loc["root:Environmental:Soil", "top_level"] == "Environmental"
    assert leaf_df.loc["root:Host-associated:Plants:Soil", "top_level"] == (
        "Host-associated"
    )
    assert result["duplicate_leaf_labels"]["Soil"]["paths"] == [
        "root:Environmental:Soil",
        "root:Host-associated:Plants:Soil",
    ]


def test_format_report_can_be_written_to_text_file(tmp_path):
    module = load_module()
    result = module.analyse_biomes(
        [
            (
                "MGF0001",
                (
                    "ids,labels,parents,counts\n"
                    "root,root,,1\n"
                    "root:Environmental,Environmental,root,1\n"
                    "root:Environmental:Soil,Soil,root:Environmental,1\n"
                ),
            )
        ]
    )
    report_path = tmp_path / "biome_report.txt"

    text = module.format_report(result, Path("figure.png"), figure_written=True)
    module.write_text_report(report_path, text)

    assert "Figure saved -> figure.png" in report_path.read_text()
    assert "Exclusive leaf paths" in report_path.read_text()


def test_format_report_includes_narrowest_families():
    module = load_module()
    result = module.analyse_biomes(
        [
            (
                "MGF0001",
                (
                    "ids,labels,parents,counts\n"
                    "root,root,,1\n"
                    "root:Environmental,Environmental,root,1\n"
                    "root:Environmental:Soil,Soil,root:Environmental,1\n"
                ),
            ),
            (
                "MGF0002",
                (
                    "ids,labels,parents,counts\n"
                    "root,root,,2\n"
                    "root:Environmental,Environmental,root,2\n"
                    "root:Environmental:Soil,Soil,root:Environmental,1\n"
                    "root:Environmental:Lake,Lake,root:Environmental,1\n"
                ),
            ),
        ]
    )

    text = module.format_report(result, Path("figure.png"), figure_written=False)

    narrowest = text.index("=== Narrowest families (least leaf biomes) ===")
    assert text.index("  MGF0001", narrowest) < text.index("  MGF0002", narrowest)


def test_parse_args_requires_db_path():
    module = load_module()

    try:
        module.parse_args([])
    except SystemExit as error:
        assert error.code == 2
    else:
        raise AssertionError("parse_args accepted a missing db_path")


def test_default_output_paths_are_under_output_directory():
    module = load_module()

    args = module.parse_args(["families.sqlite3"])
    output_figure, output_report = module.output_paths(args)

    assert args.db_path == Path("families.sqlite3")
    assert output_figure == Path("output/mgnifams_biome_distribution.png")
    assert output_report == Path("output/mgnifams_biome_distribution.txt")
