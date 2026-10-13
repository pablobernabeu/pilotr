"""How `Dataset.to_csv` writes a data set: in UTF-8, with a field quoted when it needs to be, and
with LF line endings on every platform.

The method wrote the locale encoding and joined fields with bare commas. On Windows, where the
locale encoding is usually cp1252, a Polish label raised `UnicodeEncodeError` and a Spanish one
was written as cp1252 bytes. On every platform, a comma inside a label or a column name shifted
the fields that followed it. CI runs this suite with
``-X warn_default_encoding -W error::EncodingWarning`` from Python 3.10, which fails any file
opened without an encoding, on Linux as well as on Windows.
"""
import csv, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pilotr import simulate


def _two_group(levels, name="group"):
    return {
        "name": "csv", "seed": 5,
        "units": {"subject": {"n": 6}},
        "factors": [{"name": name, "levels": levels,
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 0, "coefficients": {"effect": 0.5}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }


def _read_back(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.reader(f))


def test_non_ascii_labels_are_written_in_utf8(tmp_path):
    for levels in (["fácil", "difícil"], ["Łatwy", "Trudny"]):
        data = simulate(_two_group(levels))
        path = tmp_path / "data.csv"
        data.to_csv(path)
        rows = _read_back(path)
        assert rows[0] == data.columns
        assert [r[1] for r in rows[1:]] == data.column("group")


def test_a_comma_or_a_quote_in_a_label_or_a_column_name_stays_in_its_field(tmp_path):
    data = simulate(_two_group(["low, fast", 'high "slow"'], name="speed, accuracy"))
    path = tmp_path / "data.csv"
    data.to_csv(path)
    rows = _read_back(path)
    assert {len(r) for r in rows} == {len(data.columns)}
    assert rows[0] == data.columns
    assert [r[1] for r in rows[1:]] == data.column("speed, accuracy")


def test_numbers_read_back_exactly_and_lines_end_in_lf(tmp_path):
    data = simulate(_two_group(["a", "b"]))
    path = tmp_path / "data.csv"
    data.to_csv(path)
    raw = path.read_bytes()
    # A plain data set comes out as it always did: no quotes, LF endings, a header line first.
    assert raw.startswith(b"subject,group,y\n")
    assert b"\r" not in raw and b'"' not in raw
    rows = _read_back(path)
    assert [int(r[0]) for r in rows[1:]] == data.column("subject")
    assert [float(r[2]) for r in rows[1:]] == data.column("y")
