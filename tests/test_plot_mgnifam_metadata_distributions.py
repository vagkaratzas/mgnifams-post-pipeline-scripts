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


def test_split_groups_cuts_are_half_open_so_no_family_is_lost_or_double_counted(tmp_path):
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

    groups = module.split_groups(
        df, "size", "x", "size", str(tmp_path / "p"),
        small_max=1000, medium_max=50000,
        bins=(50, 2500, 200000), labels=("small", "medium", "large"),
    )
    sizes = [sorted(group[0]["id"].to_list()) for group in groups]

    assert sizes == [[1], [2, 3], [4]]
