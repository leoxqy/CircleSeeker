"""Unit tests for Cresil external tool wrapper."""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch, PropertyMock

import pytest

from circleseeker.external.cresil import Cresil
from circleseeker.exceptions import ExternalToolError


class TestCresilInit:
    """Tests for Cresil initialization."""

    def test_init_sets_tool_name(self):
        """Cresil should have correct tool_name."""
        with patch.object(Cresil, "_check_installation"):
            cresil = Cresil()
        assert cresil.tool_name == "cresil"

    def test_init_with_custom_logger(self):
        """Cresil should accept custom logger."""
        mock_logger = MagicMock()
        with patch.object(Cresil, "_check_installation"):
            cresil = Cresil(logger=mock_logger)
        assert cresil.logger == mock_logger

    def test_init_with_threads(self):
        """Cresil should accept threads parameter."""
        with patch.object(Cresil, "_check_installation"):
            cresil = Cresil(threads=8)
        assert cresil.threads == 8


class TestCresilTrim:
    """Tests for Cresil trim command."""

    @pytest.fixture
    def cresil(self):
        """Create Cresil instance with mocked installation check."""
        with patch.object(Cresil, "_check_installation"):
            return Cresil()

    @pytest.fixture
    def mock_files(self, tmp_path):
        """Create mock input files."""
        fasta_query = tmp_path / "query.fasta"
        fasta_query.write_text(">seq1\nACGT\n")
        reference_mmi = tmp_path / "ref.mmi"
        reference_mmi.write_text("mock mmi")
        output_dir = tmp_path / "output"
        return fasta_query, reference_mmi, output_dir

    def test_trim_builds_correct_command(self, cresil, mock_files, tmp_path):
        """Trim should build correct command line."""
        fasta_query, reference_mmi, output_dir = mock_files

        # Mock the output file creation
        trim_output = output_dir / "trim.txt"

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "")
            output_dir.mkdir(parents=True, exist_ok=True)
            trim_output.write_text("mock trim output")

            result = cresil.trim(fasta_query, reference_mmi, output_dir, threads=4)

            # Verify command was called correctly
            mock_run.assert_called_once()
            cmd = mock_run.call_args[0][0]
            assert cmd[0] == "cresil"
            assert cmd[1] == "trim"
            assert "-t" in cmd
            assert "4" in cmd
            assert "-fq" in cmd
            assert "-r" in cmd
            assert "-o" in cmd

    def test_trim_creates_output_directory(self, cresil, mock_files):
        """Trim should create output directory if it doesn't exist."""
        fasta_query, reference_mmi, output_dir = mock_files

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "")
            # Create output after run is called
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "trim.txt").write_text("mock")

            cresil.trim(fasta_query, reference_mmi, output_dir)

        assert output_dir.exists()

    def test_trim_raises_on_missing_output(self, cresil, mock_files):
        """Trim should raise FileNotFoundError if output not created."""
        fasta_query, reference_mmi, output_dir = mock_files

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "")
            output_dir.mkdir(parents=True, exist_ok=True)
            # Don't create the output file

            with pytest.raises(FileNotFoundError, match="Expected trim output not found"):
                cresil.trim(fasta_query, reference_mmi, output_dir)

    def test_trim_logs_warning_on_stderr_error(self, cresil, mock_files):
        """Trim should log warning when stderr contains ERROR."""
        fasta_query, reference_mmi, output_dir = mock_files

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "ERROR: something went wrong")
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "trim.txt").write_text("mock")

            with patch.object(cresil.logger, "warning") as mock_warning:
                cresil.trim(fasta_query, reference_mmi, output_dir)
                mock_warning.assert_called_once()


class TestCresilIdentify:
    """Tests for Cresil identify command."""

    @pytest.fixture
    def cresil(self):
        """Create Cresil instance with mocked installation check."""
        with patch.object(Cresil, "_check_installation"):
            return Cresil()

    @pytest.fixture
    def mock_files(self, tmp_path):
        """Create mock input files."""
        fasta_query = tmp_path / "query.fasta"
        fasta_query.write_text(">seq1\nACGT\n")
        reference_fasta = tmp_path / "ref.fasta"
        reference_fasta.write_text(">chr1\nACGT\n")
        reference_fai = tmp_path / "ref.fasta.fai"
        reference_fai.write_text("chr1\t4\t6\t4\t5\n")
        trim_file = tmp_path / "trim.txt"
        trim_file.write_text("mock trim")
        output_dir = tmp_path / "output"
        return fasta_query, reference_fasta, reference_fai, trim_file, output_dir

    def test_identify_builds_correct_command(self, cresil, mock_files):
        """Identify should build correct command line."""
        fasta_query, reference_fasta, reference_fai, trim_file, output_dir = mock_files

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "")
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "eccDNA_final.txt").write_text("mock output")

            cresil.identify(
                fasta_query, reference_fasta, reference_fai, trim_file, output_dir,
                threads=4, split_reads=True
            )

            mock_run.assert_called_once()
            cmd = mock_run.call_args[0][0]
            assert cmd[0] == "cresil"
            assert cmd[1] == "identify"
            assert "-t" in cmd
            assert "-fa" in cmd
            assert "-fai" in cmd
            assert "-fq" in cmd
            assert "-trim" in cmd
            assert "-s" in cmd  # split_reads flag

    def test_identify_without_split_reads(self, cresil, mock_files):
        """Identify should not include -s flag when split_reads=False."""
        fasta_query, reference_fasta, reference_fai, trim_file, output_dir = mock_files

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "")
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "eccDNA_final.txt").write_text("mock output")

            cresil.identify(
                fasta_query, reference_fasta, reference_fai, trim_file, output_dir,
                split_reads=False
            )

            cmd = mock_run.call_args[0][0]
            assert "-s" not in cmd

    def test_identify_raises_on_missing_output(self, cresil, mock_files):
        """Identify should raise FileNotFoundError if output not created."""
        fasta_query, reference_fasta, reference_fai, trim_file, output_dir = mock_files

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "")
            output_dir.mkdir(parents=True, exist_ok=True)

            with pytest.raises(FileNotFoundError, match="Expected eccDNA output not found"):
                cresil.identify(
                    fasta_query, reference_fasta, reference_fai, trim_file, output_dir
                )


class TestCresilFullPipeline:
    """Tests for Cresil full pipeline."""

    @pytest.fixture
    def cresil(self):
        """Create Cresil instance with mocked installation check."""
        with patch.object(Cresil, "_check_installation"):
            return Cresil()

    @pytest.fixture
    def mock_files(self, tmp_path):
        """Create mock input files."""
        fasta_query = tmp_path / "query.fasta"
        fasta_query.write_text(">seq1\nACGT\n")
        reference_fasta = tmp_path / "ref.fasta"
        reference_fasta.write_text(">chr1\nACGT\n")
        reference_fai = tmp_path / "ref.fasta.fai"
        reference_fai.write_text("chr1\t4\t6\t4\t5\n")
        reference_mmi = tmp_path / "ref.mmi"
        reference_mmi.write_text("mock mmi")
        output_dir = tmp_path / "output"
        return fasta_query, reference_fasta, reference_mmi, output_dir

    def test_full_pipeline_calls_trim_and_identify(self, cresil, mock_files):
        """Full pipeline should call both trim and identify."""
        fasta_query, reference_fasta, reference_mmi, output_dir = mock_files

        # Create the .fai file
        reference_fai = reference_fasta.with_suffix(reference_fasta.suffix + ".fai")
        reference_fai.write_text("chr1\t4\t6\t4\t5\n")

        with patch.object(cresil, "trim") as mock_trim, \
             patch.object(cresil, "identify") as mock_identify:

            trim_output = output_dir / "trim.txt"
            eccDNA_output = output_dir / "eccDNA_final.txt"

            mock_trim.return_value = trim_output
            mock_identify.return_value = eccDNA_output

            result = cresil.run_full_pipeline(
                fasta_query, reference_fasta, reference_mmi, output_dir
            )

            mock_trim.assert_called_once()
            mock_identify.assert_called_once()
            assert result == eccDNA_output

    def test_full_pipeline_creates_fai_if_missing(self, cresil, tmp_path):
        """Full pipeline should create .fai index if missing."""
        # Create fresh files without .fai
        fasta_query = tmp_path / "query.fasta"
        fasta_query.write_text(">seq1\nACGT\n")
        reference_fasta = tmp_path / "ref.fasta"
        reference_fasta.write_text(">chr1\nACGT\n")
        reference_mmi = tmp_path / "ref.mmi"
        reference_mmi.write_text("mock mmi")
        output_dir = tmp_path / "output"
        # Explicitly don't create .fai file

        with patch.object(cresil, "trim") as mock_trim, \
             patch.object(cresil, "identify") as mock_identify, \
             patch("circleseeker.external.cresil.Samtools") as mock_samtools_class:

            trim_output = output_dir / "trim.txt"
            eccDNA_output = output_dir / "eccDNA_final.txt"

            mock_trim.return_value = trim_output
            mock_identify.return_value = eccDNA_output

            mock_samtools = MagicMock()
            mock_samtools_class.return_value = mock_samtools

            # Make faidx create the .fai file
            def create_fai(ref_path):
                fai_path = Path(ref_path).with_suffix(Path(ref_path).suffix + ".fai")
                fai_path.write_text("chr1\t4\t6\t4\t5\n")

            mock_samtools.faidx.side_effect = create_fai

            cresil.run_full_pipeline(
                fasta_query, reference_fasta, reference_mmi, output_dir
            )

            mock_samtools.faidx.assert_called()


class TestCresilErrorHandling:
    """Tests for Cresil error handling."""

    @pytest.fixture
    def cresil(self):
        """Create Cresil instance with mocked installation check."""
        with patch.object(Cresil, "_check_installation"):
            return Cresil()

    def test_trim_handles_run_error(self, cresil, tmp_path):
        """Trim should propagate errors from run()."""
        fasta_query = tmp_path / "query.fasta"
        fasta_query.write_text(">seq1\nACGT\n")
        reference_mmi = tmp_path / "ref.mmi"
        reference_mmi.write_text("mock")
        output_dir = tmp_path / "output"

        with patch.object(cresil, "run") as mock_run:
            mock_run.side_effect = ExternalToolError("cresil failed")

            with pytest.raises(ExternalToolError):
                cresil.trim(fasta_query, reference_mmi, output_dir)

    def test_identify_handles_run_error(self, cresil, tmp_path):
        """Identify should propagate errors from run()."""
        fasta_query = tmp_path / "query.fasta"
        fasta_query.write_text(">seq1\nACGT\n")
        reference_fasta = tmp_path / "ref.fasta"
        reference_fasta.write_text(">chr1\nACGT\n")
        reference_fai = tmp_path / "ref.fasta.fai"
        reference_fai.write_text("chr1\t4\t6\t4\t5\n")
        trim_file = tmp_path / "trim.txt"
        trim_file.write_text("mock")
        output_dir = tmp_path / "output"

        with patch.object(cresil, "run") as mock_run:
            mock_run.side_effect = ExternalToolError("cresil failed")

            with pytest.raises(ExternalToolError):
                cresil.identify(
                    fasta_query, reference_fasta, reference_fai, trim_file, output_dir
                )


class TestCresilPathHandling:
    """Tests for Cresil path handling."""

    @pytest.fixture
    def cresil(self):
        """Create Cresil instance with mocked installation check."""
        with patch.object(Cresil, "_check_installation"):
            return Cresil()

    def test_trim_resolves_paths(self, cresil, tmp_path):
        """Trim should resolve all paths to absolute."""
        fasta_query = tmp_path / "query.fasta"
        fasta_query.write_text(">seq1\nACGT\n")
        reference_mmi = tmp_path / "ref.mmi"
        reference_mmi.write_text("mock")
        output_dir = tmp_path / "output"

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "")
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "trim.txt").write_text("mock")

            cresil.trim(fasta_query, reference_mmi, output_dir)

            cmd = mock_run.call_args[0][0]
            # All paths should be absolute
            for i, arg in enumerate(cmd):
                if arg in ["-fq", "-r", "-o"]:
                    path_arg = cmd[i + 1]
                    assert Path(path_arg).is_absolute()

    def test_identify_uses_cwd_for_output(self, cresil, tmp_path):
        """Identify should run from output_dir to control output location."""
        fasta_query = tmp_path / "query.fasta"
        fasta_query.write_text(">seq1\nACGT\n")
        reference_fasta = tmp_path / "ref.fasta"
        reference_fasta.write_text(">chr1\nACGT\n")
        reference_fai = tmp_path / "ref.fasta.fai"
        reference_fai.write_text("chr1\t4\t6\t4\t5\n")
        trim_file = tmp_path / "trim.txt"
        trim_file.write_text("mock")
        output_dir = tmp_path / "output"

        with patch.object(cresil, "run") as mock_run:
            mock_run.return_value = ("", "")
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "eccDNA_final.txt").write_text("mock")

            cresil.identify(
                fasta_query, reference_fasta, reference_fai, trim_file, output_dir
            )

            # Verify cwd was set
            call_kwargs = mock_run.call_args[1]
            assert "cwd" in call_kwargs
            assert Path(call_kwargs["cwd"]).resolve() == output_dir.resolve()
