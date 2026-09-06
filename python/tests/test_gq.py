from __future__ import annotations

import json
import subprocess
import unittest
from unittest.mock import patch

from gravlax import Client


class GqClientTests(unittest.TestCase):
    @patch("gravlax.client.subprocess.run")
    def test_validate_keeps_paths_as_data_and_scan_permission_explicit(self, run):
        run.return_value = subprocess.CompletedProcess([], 0, json.dumps({"valid": True}), "")
        self.assertTrue(Client().gq_validate(
            "--unusual query.gq", bindings={"pbmc": "archive with spaces.aie"}
        )["valid"])
        argv = run.call_args.args[0]
        self.assertIn("--bind=pbmc=archive with spaces.aie", argv)
        self.assertNotIn("--allow-full-scan", argv)
        self.assertNotIn("--parallel-decode", argv)
        self.assertEqual(argv[-2:], ("--", "--unusual query.gq"))
        self.assertIs(run.call_args.kwargs["shell"], False)

    def test_invalid_options_are_rejected_before_subprocess(self):
        for options in ({"max_records": 0}, {"bindings": {"x=y": "a.aie"}}, {"allow_full_scan": "yes"}, {"parallel_decode": "yes"}):
            with self.assertRaises((ValueError, TypeError)):
                Client().gq_run("query.gq", **options)


if __name__ == "__main__":
    unittest.main()
