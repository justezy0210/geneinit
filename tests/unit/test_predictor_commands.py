from pathlib import Path
import unittest

from geneinit.predictor.commands import (
    augustus_etraining_command,
    augustus_new_species_command,
    augustus_optimize_command,
    augustus_prediction_command,
    augustus_training_genbank_command,
    genemark_etp_command,
)


class PredictorCommandTests(unittest.TestCase):
    def test_genemark_is_taxon_neutral_and_preserves_softmasking(self) -> None:
        command = genemark_etp_command(
            config=Path("config.yaml"),
            workdir=Path("work"),
            bam_dir=Path("bam"),
            proc="protein_evidence",
            threads=36,
        )
        self.assertEqual(command[:2], ("/bin/bash", "-c"))
        self.assertIn('PATH="/usr/bin:/bin:${PATH}"', command[2])
        self.assertEqual(command[3:5], ("geneinit-gmetp", "/opt/ETP/bin/gmetp.pl"))
        self.assertIn("--softmask", command)
        self.assertNotIn("--fungus", command)
        self.assertNotIn("--force", command)

    def test_augustus_training_is_direct(self) -> None:
        new_species = augustus_new_species_command(
            species="ASM_etp", augustus_config=Path("config")
        )
        convert = augustus_training_genbank_command(
            training_gtf=Path("training.gtf"),
            genome=Path("genome.fa"),
            flank=1000,
            output=Path("training.gb"),
        )
        etraining = augustus_etraining_command(
            species="ASM_etp",
            training_genbank=Path("training.gb"),
            augustus_config=Path("config"),
        )
        optimize = augustus_optimize_command(
            species="ASM_etp",
            training_genbank=Path("training.gb"),
            augustus_config=Path("config"),
            threads=8,
            rounds=5,
            kfold=8,
            onlytrain=Path("onlytrain.gb"),
        )
        combined = new_species + convert + etraining + optimize
        self.assertTrue(any(value.endswith("new_species.pl") for value in combined))
        self.assertTrue(any(value.endswith("gff2gbSmallDNA.pl") for value in combined))
        self.assertTrue(any(value.endswith("etraining") for value in combined))
        self.assertTrue(any(value.endswith("optimize_augustus.pl") for value in combined))
        self.assertNotIn("--softmasked", convert)
        self.assertIn("--stopCodonExcludedFromCDS=false", etraining)
        self.assertIn("--onlytrain=onlytrain.gb", optimize)
        self.assertFalse(any("BRAKER" in value or "fungus" in value for value in combined))

    def test_augustus_prediction_uses_hints_directly(self) -> None:
        command = augustus_prediction_command(
            species="ASM_etp",
            genome=Path("genome.fa"),
            hints=Path("hints.gff"),
            extrinsic_config=Path("geneinit_etp.cfg"),
            augustus_config=Path("config"),
            output=Path("prediction.gff3"),
        )
        self.assertIn("--hintsfile=hints.gff", command)
        self.assertIn("--extrinsicCfgFile=geneinit_etp.cfg", command)
        self.assertIn("--softmasking=1", command)
        self.assertIn("--UTR=off", command)
        self.assertFalse(any("BRAKER" in value or "fungus" in value for value in command))


if __name__ == "__main__":
    unittest.main()
