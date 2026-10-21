from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from geneinit.contracts.deployment import (
    DeploymentError,
    PredictorDeployment,
    apptainer_command,
    bind_arguments,
    inspect_container,
    load_predictor_deployment,
    validate_predictor_versions,
)


class DeploymentContractTests(unittest.TestCase):
    def test_sif_checksum_and_license_are_required(self) -> None:
        with tempfile.TemporaryDirectory() as raw_tmp:
            sif = Path(raw_tmp) / "predictor.sif"
            sif.write_bytes(b"fixture")
            digest = hashlib.sha256(b"fixture").hexdigest()
            config = {
                "deployment": {
                    "predictor_sif": str(sif),
                    "predictor_sif_sha256": digest,
                    "genemark_license_acknowledged": True,
                    "apptainer_executable": sys.executable,
                }
            }
            result = load_predictor_deployment(config)
            self.assertEqual(result.mode, "apptainer")
            self.assertEqual(result.sif_sha256, digest)

            config["deployment"]["predictor_sif_sha256"] = "0" * 64
            with self.assertRaisesRegex(DeploymentError, "mismatch"):
                load_predictor_deployment(config)

    def test_missing_license_acknowledgement_fails(self) -> None:
        with self.assertRaisesRegex(DeploymentError, "license"):
            load_predictor_deployment({"deployment": {}})

    def test_in_container_mode_does_not_require_apptainer(self) -> None:
        with tempfile.TemporaryDirectory() as raw_tmp:
            root = Path(raw_tmp)
            sif = root / "allinone.sif"
            sif.write_bytes(b"fixture")
            config = {
                "deployment": {
                    "mode": "in_container",
                    "predictor_sif": str(sif),
                    "predictor_sif_sha256": hashlib.sha256(b"fixture").hexdigest(),
                    "genemark_license_acknowledged": True,
                }
            }
            result = load_predictor_deployment(config)
            self.assertEqual(result.mode, "in_container")
            self.assertIsNone(result.apptainer)
            command = apptainer_command(result, ("tool", "--flag"), (), pwd=root)
            self.assertEqual(command[:2], ["/usr/bin/env", f"--chdir={root.resolve()}"])
            self.assertEqual(command[-2:], ["tool", "--flag"])
            self.assertIn("OPENBLAS_NUM_THREADS=1", command)
            self.assertIn("OMP_NUM_THREADS=1", command)

    def test_bind_roots_are_deduplicated(self) -> None:
        paths = [Path("/tmp/run/a.txt"), Path("/tmp/run/sub/b.txt")]
        self.assertEqual(bind_arguments(paths), ["--bind", "/tmp/run:/tmp/run"])

    def test_existing_input_file_is_bound_as_a_file(self) -> None:
        with tempfile.TemporaryDirectory() as raw_tmp:
            source = Path(raw_tmp) / "private.sif"
            source.write_bytes(b"fixture")
            self.assertEqual(
                bind_arguments([source]),
                ["--bind", f"{source.resolve()}:{source.resolve()}"],
            )

    def test_container_inspection_checks_system_perl_runtime(self) -> None:
        deployment = PredictorDeployment(
            mode="in_container",
            sif=Path("/tmp/predictor.sif"),
            sif_sha256="0" * 64,
            apptainer=None,
        )
        completed = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout='AUGUSTUS (3.5.0)\nmy $version = "1.02";\n',
            stderr="",
        )
        with patch(
            "geneinit.contracts.deployment.subprocess.run", return_value=completed
        ) as run:
            inspect_container(deployment)

        command = run.call_args.args[0]
        self.assertEqual(command[-3:-1], ["/bin/bash", "-lc"])
        self.assertIn("/usr/bin/perl -MHash::Merge", command[-1])
        self.assertIn("-MYAML -MYAML::XS", command[-1])
        self.assertIn("find /opt/ETP/bin", command[-1])

    def test_predictor_versions_are_exact(self) -> None:
        text = (
            "AUGUSTUS (3.5.0)\n"
            'my $version = "1.02";\n'
        )
        observed = validate_predictor_versions(text)
        self.assertEqual(observed["GeneMark-ETP"], "1.02")
        with self.assertRaisesRegex(DeploymentError, "unsupported AUGUSTUS"):
            validate_predictor_versions(text.replace("3.5.0", "3.6.0"))


if __name__ == "__main__":
    unittest.main()
