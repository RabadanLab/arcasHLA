"""
Direct code testing of the reference customization script.
"""

import os.path

import customize


class TestMain:
    def test_basic(self, repo_root, tmp_path):
        genotype_result_path = os.path.join(
            repo_root, "test/expected_output/test.genotype.json"
        )

        custom_ref_out_path = os.path.join(str(tmp_path), "custom_reference")

        main_args = [
            "--genotype",
            genotype_result_path,
            "-o",
            custom_ref_out_path,
        ]

        # Should run without error.
        customize.main(main_args)
