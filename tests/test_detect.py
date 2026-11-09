"""Unit tests for the pure functions in debreak_detect.py.

These are characterization tests: they pin the current behaviour.
Places where the current behaviour looks like a bug are marked
`SUSPECTED BUG` so the refactor can decide consciously whether to fix them.
"""
import pytest

import debreak_detect as dd


# ---------------------------------------------------------------- cigardeletion

def _cd(cigar, flag=0, pos=1000, min_size=48, max_size=400_000_000):
    return dd.cigardeletion(flag, "chr1", pos, cigar, min_size, max_size)


def test_cigardeletion_single_deletion():
    svs, reflen, clips = _cd("100M200D100M")
    assert svs == [["chr1", 1100, 200, "D-cigar"]]
    assert reflen == 400
    assert clips == [0, 200, 0]


def test_cigardeletion_single_insertion_records_read_offset():
    svs, reflen, clips = _cd("100M60I100M")
    assert svs == [["chr1", 1100, 60, "I-cigar", 100]]
    assert reflen == 200
    assert clips == [0, 260, 0]


def test_cigardeletion_read_offset_includes_left_soft_clip():
    svs, reflen, clips = _cd("50S100M60I100M30S")
    assert svs == [["chr1", 1100, 60, "I-cigar", 150]]
    assert clips == [50, 260, 30]


def test_cigardeletion_hard_clips_counted_like_soft_clips():
    _, _, clips = _cd("40H100M25H")
    assert clips == [40, 100, 25]


@pytest.mark.parametrize("op", ["M", "=", "X"])
def test_cigardeletion_match_ops_consume_ref_and_read(op):
    _, reflen, clips = _cd(f"100{op}")
    assert reflen == 100
    assert clips == [0, 100, 0]


def test_cigardeletion_respects_size_bounds():
    assert _cd("100M47D100M")[0] == []  # below default min_size=48
    assert _cd("100M48D100M")[0] == [["chr1", 1100, 48, "D-cigar"]]
    assert _cd("100M300D100M", max_size=299)[0] == []
    assert _cd("100M300I100M", max_size=299)[0] == []


def test_cigardeletion_small_events_still_advance_coordinates():
    svs, reflen, _ = _cd("100M10D100M200D100M")
    assert svs == [["chr1", 1210, 200, "D-cigar"]]
    assert reflen == 510


@pytest.mark.parametrize("flag", [0, 16])
def test_cigardeletion_primary_flags_report_svs(flag):
    assert _cd("100M200D100M", flag=flag)[0] != []


@pytest.mark.parametrize("flag", [256, 272, 2048, 2064])
def test_cigardeletion_non_primary_flags_report_nothing(flag):
    svs, reflen, clips = _cd("100M200D100M", flag=flag)
    assert svs == []
    assert reflen == 400  # coordinates are still computed
    assert clips == [0, 200, 0]


def test_cigardeletion_flag_may_be_string():
    assert _cd("100M200D100M", flag="16")[0] == [["chr1", 1100, 200, "D-cigar"]]


def test_cigardeletion_close_deletions_are_merged_by_summing_sizes():
    # gap between deletions = 300 <= 500 -> merged.
    # NOTE: merged size is the SUM of sizes (200), not the span (500).
    svs, _, _ = _cd("100M100D300M100D100M")
    assert svs == [["chr1", 1100, 200, "D-cigar"]]


def test_cigardeletion_distant_deletions_are_not_merged():
    svs, _, _ = _cd("100M100D600M100D100M")
    assert svs == [["chr1", 1100, 100, "D-cigar"], ["chr1", 1800, 100, "D-cigar"]]


def test_cigardeletion_three_close_deletions_collapse_to_one():
    svs, _, _ = _cd("100M60D100M60D100M60D100M")
    assert svs == [["chr1", 1100, 180, "D-cigar"]]


def test_cigardeletion_large_insertions_merge_within_600():
    svs, _, _ = _cd("100M600I550M600I100M")
    assert svs == [["chr1", 1100, 1200, "I-cigar", 100]]


def test_cigardeletion_large_insertions_beyond_600_not_merged():
    svs, _, _ = _cd("100M600I650M600I100M")
    assert [s[1] for s in svs] == [1100, 1750]


def test_cigardeletion_small_insertion_window_is_600_not_200():
    # SUSPECTED BUG: the code intends window=200 for insertions <100bp
    # (`window=200 if max(l1,l2)<100 else 400`), but the next line
    # `window=400 if window==400 and ... else 600` overwrites 200 with 600.
    # Two 60bp insertions 500bp apart are therefore merged.
    svs, _, _ = _cd("100M60I500M60I100M")
    assert svs == [["chr1", 1100, 120, "I-cigar", 100]]


def test_cigardeletion_mid_insertion_window_is_400():
    # max size in [100, 500) -> window 400
    merged, _, _ = _cd("100M200I400M200I100M")
    assert merged == [["chr1", 1100, 400, "I-cigar", 100]]
    split, _, _ = _cd("100M200I450M200I100M")
    assert len(split) == 2


def test_cigardeletion_deletions_listed_before_insertions():
    svs, _, _ = _cd("100M60I1000M200D100M")
    assert [s[3] for s in svs] == ["D-cigar", "I-cigar"]


# ---------------------------------------------------------------- simplifycigar

@pytest.mark.parametrize(
    "cigar, expected",
    [
        ("100M", [0, 100, 0]),
        ("50S100M", [50, 100, 0]),
        ("100M30S", [0, 100, 30]),
        ("50S100M10I20D30M40H", [50, 140, 40]),  # I counts, D does not
        ("10H5S100M", [5, 100, 0]),  # SUSPECTED QUIRK: last leading clip wins
    ],
)
def test_simplifycigar(cigar, expected):
    assert dd.simplifycigar(cigar) == expected


def test_simplifycigar_matches_cigardeletion_clip_info():
    cigar = "50S100M60I100M200D10M30S"
    assert dd.simplifycigar(cigar) == _cd(cigar)[2]


# ---------------------------------------------------------------- segmentdeletion
# readinfo layout: [name, flag, chrom, ref_start, ref_end,
#                   [left_clip, aligned_read_len, right_clip], mapq, read_seq]

SEQ = "ACGT" * 1000  # 4000bp read


def _seg(flag, chrom, start, end, clips, mapq=60, seq=""):
    return ["r1", flag, chrom, start, end, clips, mapq, seq]


PRIMARY = _seg(0, "chr1", 1000, 2000, [0, 1000, 1000], seq=SEQ)


def _sd(*supps, primary=PRIMARY, min_size=48, max_size=400_000_000):
    return dd.segmentdeletion([primary, *supps], min_size, max_size)


def test_segmentdeletion_deletion():
    supp = _seg(2048, "chr1", 2500, 3500, [1000, 1000, 0])
    assert _sd(supp) == ["chr1\t2000\t500\tD-segment\tr1_seg1\t2048\t60"]


def test_segmentdeletion_deletion_is_symmetric_in_segment_order():
    supp = _seg(2048, "chr1", 2500, 3500, [1000, 1000, 0])
    assert dd.segmentdeletion([supp, PRIMARY], 48, 400_000_000) == _sd(supp)


def test_segmentdeletion_insertion_includes_inserted_sequence():
    primary = _seg(0, "chr1", 1000, 2000, [0, 1000, 1500], seq=SEQ[:2500])
    supp = _seg(2048, "chr1", 2000, 3000, [1500, 1000, 0])
    assert _sd(supp, primary=primary) == [
        "chr1\t2000\t500\tI-segment\tr1_seg1\t2048\t60\t" + SEQ[1000:1500]
    ]


def test_segmentdeletion_inversion():
    supp = _seg(2064, "chr1", 2500, 3500, [0, 1000, 1000])
    assert _sd(supp) == ["chr1\t2000\t1500\tINV-segment\tr1_seg1\t2064\t60"]


def test_segmentdeletion_translocation():
    supp = _seg(2048, "chr2", 5000, 6000, [1000, 1000, 0])
    assert _sd(supp) == ["chr1\t2000\tchr2\t5000\tTRA-segment\tr1_seg1\t2048\t60"]


def test_segmentdeletion_translocation_orders_chromosomes():
    primary = _seg(0, "chr2", 1000, 2000, [0, 1000, 1000], seq=SEQ)
    supp = _seg(2048, "chr1", 5000, 6000, [1000, 1000, 0])
    assert _sd(supp, primary=primary) == [
        "chr1\t5000\tchr2\t2000\tTRA-segment\tr1_seg1\t2048\t60"
    ]


def test_segmentdeletion_mapq_is_mean_of_segments():
    supp = _seg(2048, "chr1", 2500, 3500, [1000, 1000, 0], mapq=21)
    assert _sd(supp)[0].endswith("\t40")  # (60 + 21) // 2


def test_segmentdeletion_ids_increment_per_call():
    supp1 = _seg(2048, "chr1", 2500, 3500, [1000, 1000, 0])
    supp2 = _seg(2048, "chr2", 5000, 6000, [1000, 1000, 0])
    calls = _sd(supp1, supp2)
    assert [c.split("\t")[-3] for c in calls] == ["r1_seg1", "r1_seg2"]


def test_segmentdeletion_respects_size_bounds():
    supp = _seg(2048, "chr1", 2500, 3500, [1000, 1000, 0])
    assert _sd(supp, min_size=501) == []
    assert _sd(supp, max_size=499) == []


def test_segmentdeletion_ignores_short_segments():
    supp = _seg(2048, "chr1", 2500, 2799, [1000, 299, 0])  # aligned < 300bp
    assert _sd(supp) == []


def test_segmentdeletion_needs_a_primary_alignment():
    supp1 = _seg(2048, "chr1", 1000, 2000, [0, 1000, 1000])
    supp2 = _seg(2048, "chr1", 2500, 3500, [1000, 1000, 0])
    assert dd.segmentdeletion([supp1, supp2], 48, 400_000_000) == []


def test_segmentdeletion_single_segment_gives_nothing():
    assert dd.segmentdeletion([PRIMARY], 48, 400_000_000) == []


def test_segmentdeletion_same_start_is_skipped():
    supp = _seg(2048, "chr1", 1000, 2000, [1000, 1000, 0])
    assert _sd(supp) == []


# ---------------------------------------------------------------- segmentdeletion_tumor

def test_segmentdeletion_tumor_matches_normal_for_deletion():
    supp = _seg(2048, "chr1", 2500, 3500, [1000, 1000, 0])
    assert dd.segmentdeletion_tumor([PRIMARY, supp], 48, 400_000_000) == _sd(supp)


def test_segmentdeletion_tumor_insertion_omits_sequence():
    # Difference vs segmentdeletion: no inserted sequence column.
    primary = _seg(0, "chr1", 1000, 2000, [0, 1000, 1500], seq=SEQ[:2500])
    supp = _seg(2048, "chr1", 2000, 3000, [1500, 1000, 0])
    assert dd.segmentdeletion_tumor([primary, supp], 48, 400_000_000) == [
        "chr1\t2000\t500\tI-segment\tr1_seg1\t2048\t60"
    ]


def test_segmentdeletion_tumor_compares_all_pairs_without_primary():
    s1 = _seg(2048, "chr1", 1000, 2000, [0, 1000, 1000])
    s2 = _seg(2048, "chr1", 2500, 3500, [1000, 1000, 0])
    assert dd.segmentdeletion_tumor([s1, s2], 48, 400_000_000) == [
        "chr1\t2000\t500\tD-segment\tr1_seg1\t2048\t60"
    ]


# ---------------------------------------------------------------- detect_sortbam

def test_detect_sortbam_writes_raw_calls_and_depth(tmp_path):
    from conftest import REPO
    import os

    (tmp_path / "sv_raw_calls").mkdir()
    (tmp_path / "map_depth").mkdir()
    bam = os.path.join(REPO, "testdata", "test_read.bam")
    assert dd.detect_sortbam(bam, str(tmp_path) + "/", 48, 400_000_000,
                             "chr1", ["chr1"], False, False, False)

    raw = (tmp_path / "sv_raw_calls" / "test_read.bam-chr1.debreak.temp").read_text()
    golden = open(os.path.join(REPO, "tests", "golden", "default", "sv_raw_calls.txt")).read()
    assert sorted(raw.splitlines()) == golden.splitlines()
    assert (tmp_path / "map_depth" / "maplength_chr1").read_text() == "22577500\n"
