import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "bin"
    / "plot_mgnifam_metadata_distributions.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "plot_mgnifam_metadata_distributions", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_inputs(tmp_path, rows, novel_ids):
    metadata = tmp_path / "mgnifam.csv"
    header = "id,full_size,protein_rep,rep_region,rep_length,converged,plddt,ptm\n"
    body = "".join(
        f"{i},{size},999,1-10,{length},True,{plddt},0.5\n"
        for i, size, length, plddt in rows
    )
    metadata.write_text(header + body)

    novel = tmp_path / "novel_ids.txt"
    novel.write_text("".join(f"{i}\n" for i in novel_ids))
    return metadata, novel


def test_annotated_flag_is_the_complement_of_the_novel_id_list(tmp_path):
    module = load_module()
    metadata, novel = write_inputs(
        tmp_path,
        rows=[(1, 30, 100, 50.0), (2, 5000, 400, 80.0), (3, 200000, 1500, 90.0)],
        novel_ids=[2],
    )

    df = module.load_metadata(str(metadata), str(novel))

    assert dict(zip(df["id"], df["annotated"])) == {1: True, 2: False, 3: True}


def test_split_panels_cuts_are_half_open_so_no_family_is_lost_or_double_counted(tmp_path):
    module = load_module()
    metadata, novel = write_inputs(
        tmp_path,
        # 999/1000 and 49,999/50,000 sit exactly on the small/medium and medium/large cuts
        rows=[
            (1, 999, 100, 50.0),
            (2, 1000, 100, 50.0),
            (3, 49999, 100, 50.0),
            (4, 50000, 100, 50.0),
        ],
        novel_ids=[],
    )
    df = module.load_metadata(str(metadata), str(novel))

    groups = module.split_panels(
        df, "size", "x", "size",
        small_max=1000, medium_max=50000,
        bins=(50, 2500, 200000), labels=("small", "medium", "large"),
    )
    sizes = [sorted(group[0]["id"].to_list()) for group in groups]

    assert sizes == [[1], [2, 3], [4]]


def test_every_number_is_labelled_once_and_no_two_labels_on_a_bar_overlap(tmp_path):
    """The point of the measured layout: nothing dropped, nothing on top of anything else."""
    module = load_module()
    # a tall bar (labels fit inside) next to a sliver (everything is forced above the bar)
    metadata, novel = write_inputs(
        tmp_path,
        rows=[(i, 150, 100, 50.0) for i in range(1, 4001)]
        + [(i, 150, 100, 90.0) for i in range(4001, 4004)],
        novel_ids=list(range(3001, 4003)),
    )
    df = module.load_metadata(str(metadata), str(novel))
    tidy = module.bin_counts(df, value_col="plddt", bin_size=5.0)

    plot_w_mm, plot_h_mm = module.panel_area_mm(tidy, "x", module.PANEL_HEIGHT_MM)
    y_top = tidy["total"].max() * 4  # generous, so the layout is not the one clipping anything
    labels, y_reach = module.layout_labels(tidy, plot_w_mm, plot_h_mm, y_top)

    assert set(labels["label"]) == {"3,000", "1,000", "25.0%", "1", "2", "66.7%"}
    assert y_reach <= y_top

    units_per_mm = y_top / plot_h_mm
    for _, bar in labels.groupby("bin", observed=True):
        spans = []
        for row in bar.itertuples():
            w, h = module.text_size_mm(row.label)
            extent = (w if row.angle == 90 else h) * units_per_mm  # upright text runs its width
            spans.append((row.y - extent / 2, row.y + extent / 2) if row.va == "center"
                         else (row.y, row.y + extent))
        spans.sort()
        assert all(a[1] <= b[0] for a, b in zip(spans, spans[1:])), spans
