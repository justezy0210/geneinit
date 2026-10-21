import unittest

from geneinit.evidence.junctions import (
    cigar_junctions,
    recurrence_metrics,
    transcript_strand,
)


class FakeAlignment:
    def __init__(self, *, reverse=False, read1=False, read2=False, xs=None):
        self.is_reverse = reverse
        self.is_read1 = read1
        self.is_read2 = read2
        self._xs = xs

    def has_tag(self, name):
        return name == "XS" and self._xs is not None

    def get_tag(self, name):
        if not self.has_tag(name):
            raise KeyError(name)
        return self._xs


class JunctionEvidenceTests(unittest.TestCase):
    def test_cigar_coordinates_and_anchors(self):
        observed = cigar_junctions(
            99, [(4, 3), (0, 10), (3, 100), (0, 20), (3, 50), (0, 7)]
        )
        self.assertEqual(observed, [(110, 209, 10, 20), (230, 279, 20, 7)])

    def test_xs_has_priority(self):
        alignment = FakeAlignment(reverse=False, read1=True, xs="-")
        self.assertEqual(transcript_strand(alignment, "FR"), ("-", "XS"))

    def test_pe_and_se_protocol_fallbacks(self):
        read1_reverse = FakeAlignment(reverse=True, read1=True)
        read2_reverse = FakeAlignment(reverse=True, read2=True)
        self.assertEqual(transcript_strand(read1_reverse, "RF"), ("+", "inferred"))
        self.assertEqual(transcript_strand(read2_reverse, "RF"), ("-", "inferred"))
        self.assertEqual(transcript_strand(read1_reverse, "R"), ("+", "inferred"))
        self.assertEqual(transcript_strand(read1_reverse, "F"), ("-", "inferred"))

    def test_recurrence_is_condition_aware(self):
        result = recurrence_metrics(
            [3, 2, 0, 4], ["conditionA", "conditionA", "conditionA", "conditionB"]
        )
        self.assertEqual(result["sample_count_hq_ge1"], 3)
        self.assertEqual(result["condition_count_ge2rep_hq_ge1"], 1)
        self.assertEqual(result["condition_count_ge3rep_hq_ge1"], 0)


if __name__ == "__main__":
    unittest.main()
