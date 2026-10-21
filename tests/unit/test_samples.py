from pathlib import Path
import os
from tempfile import TemporaryDirectory
import unittest

from geneinit.contracts.samples import SampleContractError, load_samples


class SampleProtocolTests(unittest.TestCase):
    EXPECTED = {
        ("PE", "RF"): ("RF", "--rf"),
        ("PE", "FR"): ("FR", "--fr"),
        ("PE", "unstranded"): (None, None),
        ("SE", "R"): ("R", "--rf"),
        ("SE", "F"): ("F", "--fr"),
        ("SE", "unstranded"): (None, None),
    }

    def test_all_protocol_mappings(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            read1 = root / "R1.fastq"
            read2 = root / "R2.fastq"
            read1.write_text("@r\nACGT\n+\nFFFF\n", encoding="utf-8")
            read2.write_text("@r\nTGCA\n+\nFFFF\n", encoding="utf-8")
            for index, ((layout, strand), expected) in enumerate(
                self.EXPECTED.items(), 1
            ):
                samples = root / f"samples_{index}.tsv"
                second = str(read2) if layout == "PE" else ""
                samples.write_text(
                    "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
                    f"S{index}\tC\t{index}\t{layout}\t{strand}\t{read1}\t{second}\n",
                    encoding="utf-8",
                )
                record = load_samples(samples)[0]
                self.assertEqual(
                    (record.hisat2_strandness, record.stringtie_flag), expected
                )
                self.assertEqual(record.ingenannot_paired, layout == "PE")
                self.assertEqual(
                    record.ingenannot_stranded, strand in {"RF", "R"}
                )
                if strand in {"FR", "F"}:
                    self.assertIn("FALLBACK", record.ingenannot_protocol_mode)

    def test_layout_specific_strand_values_fail(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            read1 = root / "R1.fastq"
            read1.write_text("@r\nACGT\n+\nFFFF\n", encoding="utf-8")
            samples = root / "samples.tsv"
            samples.write_text(
                "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
                f"S1\tC\t1\tSE\tRF\t{read1}\t\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(SampleContractError, "invalid protocol"):
                load_samples(samples)

    def test_duplicate_fastq_paths_and_aliases_are_not_independent_samples(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            read1, read2 = root / "R1.fastq", root / "R2.fastq"
            read1.write_text("@r\nACGT\n+\nFFFF\n")
            read2.write_text("@r\nTGCA\n+\nFFFF\n")
            symlink, hardlink = root / "symlink.fastq", root / "hardlink.fastq"
            symlink.symlink_to(read1)
            os.link(read1, hardlink)
            for alias in (read1, symlink, hardlink):
                with self.subTest(alias=alias.name):
                    samples = root / "samples.tsv"
                    samples.write_text(
                        "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
                        f"S1\tC\t1\tPE\tRF\t{read1}\t{read2}\n"
                        f"S2\tC\t2\tSE\tR\t{alias}\t\n"
                    )
                    with self.assertRaisesRegex(SampleContractError, "duplicate FASTQ"):
                        load_samples(samples)

    def test_hardlinked_mates_fail(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            read1, read2 = root / "R1.fastq", root / "R2.fastq"
            read1.write_text("@r\nACGT\n+\nFFFF\n")
            os.link(read1, read2)
            samples = root / "samples.tsv"
            samples.write_text(
                "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
                f"S1\tC\t1\tPE\tRF\t{read1}\t{read2}\n"
            )
            with self.assertRaisesRegex(SampleContractError, "read1 and read2 are identical"):
                load_samples(samples)


if __name__ == "__main__":
    unittest.main()
