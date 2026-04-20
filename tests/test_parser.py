from pathlib import Path

from worker.parser import build_count_matrix, parse_sample_sheet


def test_build_count_matrix(tmp_path: Path) -> None:
    fixture_root = Path(__file__).parent / "fixtures"
    records = parse_sample_sheet(fixture_root / "sample_sheet.csv")
    out_csv = tmp_path / "counts.csv"
    n_samples, n_genes = build_count_matrix(records, fixture_root / "rsem", out_csv)

    assert n_samples == 4
    assert n_genes == 2
    text = out_csv.read_text(encoding="utf-8")
    assert "ENSG000001.1" in text
    assert "ENSG000002.1" in text
