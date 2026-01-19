"""Tests for external Minimap2 tools."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock, mock_open
import shutil

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.minimap2 import Minimap2, MiniMapConfig


class TestMiniMapConfig:
    """Test cases for MiniMapConfig."""

    def test_minimap_config_defaults(self):
        """Test MiniMapConfig default values."""
        config = MiniMapConfig()
        assert config.preset == "map-hifi"
        assert config.threads == 8
        assert config.output_format == "bam"
        assert config.allow_secondary is True
        assert config.build_index is True
        assert config.sort_bam is True
        assert config.index_bam is True
        assert config.additional_args == ""

    def test_minimap_config_custom(self):
        """Test MiniMapConfig with custom values."""
        config = MiniMapConfig(
            preset="map-ont",
            threads=4,
            output_format="sam",
            allow_secondary=False,
            build_index=False,
            sort_bam=False,
            index_bam=False,
            additional_args="-k 15"
        )
        assert config.preset == "map-ont"
        assert config.threads == 4
        assert config.output_format == "sam"
        assert config.allow_secondary is False
        assert config.build_index is False
        assert config.sort_bam is False
        assert config.index_bam is False
        assert config.additional_args == "-k 15"


class TestMinimap2:
    """Test cases for Minimap2."""

    @patch('circleseeker.external.minimap2.shutil.which')
    def test_minimap2_initialization(self, mock_which):
        """Test Minimap2 initialization."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        minimap2 = Minimap2()
        assert minimap2.tool_name == "minimap2"
        assert minimap2.threads == 8  # Default config threads
        assert minimap2.preset == "map-hifi"
        assert minimap2.has_samtools is True

    @patch('circleseeker.external.minimap2.shutil.which')
    def test_minimap2_custom_config(self, mock_which):
        """Test Minimap2 with custom config."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        config = MiniMapConfig(threads=4, preset="map-ont")
        minimap2 = Minimap2(config=config)
        assert minimap2.threads == 4
        assert minimap2.preset == "map-ont"

    @patch('circleseeker.external.minimap2.shutil.which')
    def test_minimap2_missing_tools(self, mock_which):
        """Test Minimap2 with missing tools."""
        mock_which.return_value = None  # No tools found

        with pytest.raises(Exception):  # Should raise PipelineError
            Minimap2()

    @patch('circleseeker.external.minimap2.shutil.which')
    def test_minimap2_missing_samtools(self, mock_which):
        """Test Minimap2 with missing samtools."""
        mock_which.side_effect = lambda tool: "/usr/bin/minimap2" if tool == "minimap2" else None

        minimap2 = Minimap2()
        assert minimap2.has_samtools is False
        assert minimap2.config.sort_bam is False
        assert minimap2.config.index_bam is False

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('subprocess.run')
    def test_build_index(self, mock_run, mock_which, tmp_path):
        """Test build_index method."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )
        mock_run.return_value = MagicMock(returncode=0, stderr="")

        reference_fasta = tmp_path / "reference.fasta"
        reference_fasta.write_text(">chr1\nATCGATCG\n")

        minimap2 = Minimap2()
        index_path = minimap2.build_index(reference_fasta)

        # Verify subprocess was called
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "minimap2" in call_args
        assert "-d" in call_args
        assert str(reference_fasta) in call_args

        # Check index path generation
        assert index_path.suffix == ".mmi"

    @patch('circleseeker.external.minimap2.shutil.which')
    def test_build_index_existing(self, mock_which, tmp_path):
        """Test build_index with existing index."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        reference_fasta = tmp_path / "reference.fasta"
        index_path = tmp_path / "reference.mmi"

        # Create existing index
        index_path.touch()

        minimap2 = Minimap2()
        result_path = minimap2.build_index(reference_fasta, index_path)

        assert result_path == index_path

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('subprocess.run')
    def test_build_index_failure(self, mock_run, mock_which, tmp_path):
        """Test build_index failure handling."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        from subprocess import CalledProcessError
        mock_run.side_effect = CalledProcessError(1, "minimap2", stderr="Error building index")

        reference_fasta = tmp_path / "reference.fasta"
        reference_fasta.write_text(">chr1\nATCG\n")

        minimap2 = Minimap2()

        with pytest.raises(Exception):  # Should raise PipelineError
            minimap2.build_index(reference_fasta)

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('circleseeker.external.minimap2.Minimap2._run_alignment_simple')
    def test_align_simple(self, mock_align, mock_which, tmp_path):
        """Test align method with simple alignment."""
        mock_which.side_effect = lambda tool: f"/usr/bin/{tool}" if tool == "minimap2" else None

        reference = tmp_path / "reference.fasta"
        query = tmp_path / "query.fasta"
        output = tmp_path / "output.sam"

        reference.write_text(">ref\nATCG\n")
        query.write_text(">query\nATCG\n")

        config = MiniMapConfig(output_format="sam", build_index=False, sort_bam=False, index_bam=False)
        minimap2 = Minimap2(config=config)

        result = minimap2.align(reference, query, output)

        mock_align.assert_called_once()
        assert result == output

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('circleseeker.external.minimap2.Minimap2._run_alignment_safe')
    @patch('circleseeker.external.minimap2.Minimap2._index_bam')
    @patch('circleseeker.external.minimap2.Minimap2._print_stats')
    def test_align_bam_with_processing(self, mock_stats, mock_index, mock_align, mock_which, tmp_path):
        """Test align method with BAM processing."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        reference = tmp_path / "reference.fasta"
        query = tmp_path / "query.fasta"
        output = tmp_path / "output.bam"

        reference.write_text(">ref\nATCG\n")
        query.write_text(">query\nATCG\n")

        minimap2 = Minimap2(config=MiniMapConfig(build_index=False))
        result = minimap2.align(reference, query, output)

        mock_align.assert_called_once()
        mock_index.assert_called_once_with(output)
        mock_stats.assert_called_once_with(output)
        assert result == output

    @patch('circleseeker.external.minimap2.shutil.which')
    def test_align_missing_query(self, mock_which, tmp_path):
        """Test align with missing query file."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        reference = tmp_path / "reference.fasta"
        query = tmp_path / "missing.fasta"
        output = tmp_path / "output.sam"

        reference.write_text(">ref\nATCG\n")

        minimap2 = Minimap2()

        with pytest.raises(FileNotFoundError):
            minimap2.align(reference, query, output)

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('subprocess.run')
    def test_run_alignment_safe(self, mock_run, mock_which, tmp_path):
        """Test _run_alignment_safe method."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        # Mock successful runs for both minimap2 and samtools
        def run_side_effect(cmd, **kwargs):
            # Create temp sam file when minimap2 runs
            if "minimap2" in cmd[0]:
                temp_sam = tmp_path / "tmp_sort" / "temp_alignment.sam"
                temp_sam.parent.mkdir(exist_ok=True)
                temp_sam.write_text("@HD\tVN:1.6\n")
            # Create output bam when samtools runs
            if "samtools" in cmd[0]:
                output_bam = tmp_path / "output.bam"
                output_bam.write_bytes(b"dummy bam")
            return MagicMock(returncode=0, stderr="")

        mock_run.side_effect = run_side_effect

        reference = tmp_path / "reference.fasta"
        reads = tmp_path / "reads.fasta"
        output_bam = tmp_path / "output.bam"
        temp_dir = tmp_path / "tmp_sort"

        temp_dir.mkdir(exist_ok=True)

        minimap2 = Minimap2()
        minimap2._run_alignment_safe(reference, reads, output_bam, temp_dir)

        # Verify both minimap2 and samtools were called
        assert mock_run.call_count == 2

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('subprocess.run')
    def test_run_alignment_simple(self, mock_run, mock_which, tmp_path):
        """Test _run_alignment_simple method."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )
        mock_run.return_value = MagicMock(returncode=0, stderr="")

        reference = tmp_path / "reference.fasta"
        reads = tmp_path / "reads.fasta"
        output_file = tmp_path / "output.sam"

        minimap2 = Minimap2()
        minimap2._run_alignment_simple(reference, reads, output_file)

        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "minimap2" in call_args
        assert "-ax" in call_args
        assert str(reference) in call_args

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('subprocess.run')
    def test_additional_args_preserve_quotes_simple(self, mock_run, mock_which, tmp_path):
        """Quoted additional_args should stay grouped in simple alignment."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )
        mock_run.return_value = MagicMock(returncode=0, stderr="")

        config = MiniMapConfig(
            output_format="sam",
            build_index=False,
            sort_bam=False,
            index_bam=False,
            additional_args='--foo "bar baz"',
        )
        minimap2 = Minimap2(config=config)

        reference = tmp_path / "reference.fasta"
        reads = tmp_path / "reads.fasta"
        output_file = tmp_path / "output.sam"

        minimap2._run_alignment_simple(reference, reads, output_file)

        call_args = mock_run.call_args[0][0]
        assert "--foo" in call_args
        assert "bar baz" in call_args

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('subprocess.run')
    def test_additional_args_preserve_quotes_unsorted(self, mock_run, mock_which, tmp_path):
        """Quoted additional_args should stay grouped in unsorted BAM alignment."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        temp_dir = tmp_path / "tmp_sort"
        temp_dir.mkdir(exist_ok=True)
        output_bam = tmp_path / "output.bam"
        calls = []

        def run_side_effect(cmd, **kwargs):
            calls.append(cmd)
            if cmd[0] == "minimap2":
                temp_sam = temp_dir / "temp_alignment.sam"
                temp_sam.write_text("@HD\tVN:1.6\n")
            if cmd[0] == "samtools":
                output_bam.write_bytes(b"dummy bam")
            return MagicMock(returncode=0, stderr="")

        mock_run.side_effect = run_side_effect

        config = MiniMapConfig(
            output_format="bam",
            build_index=False,
            sort_bam=False,
            index_bam=False,
            additional_args='--foo "bar baz"',
        )
        minimap2 = Minimap2(config=config)

        reference = tmp_path / "reference.fasta"
        reads = tmp_path / "reads.fasta"

        minimap2._run_alignment_unsorted(reference, reads, output_bam, temp_dir)

        minimap2_cmd = calls[0]
        assert "--foo" in minimap2_cmd
        assert "bar baz" in minimap2_cmd

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('subprocess.run')
    def test_index_bam(self, mock_run, mock_which, tmp_path):
        """Test _index_bam method."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )
        mock_run.return_value = MagicMock(returncode=0, stderr="")

        bam_file = tmp_path / "test.bam"
        bam_file.touch()

        minimap2 = Minimap2()
        minimap2._index_bam(bam_file)

        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "samtools" in call_args
        assert "index" in call_args
        assert str(bam_file) in call_args

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('subprocess.run')
    def test_print_stats(self, mock_run, mock_which, tmp_path):
        """Test _print_stats method."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        # Mock flagstat output
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="100 reads mapped\n90 properly paired\n",
            stderr=""
        )

        bam_file = tmp_path / "test.bam"
        bam_file.write_bytes(b"dummy bam content")

        minimap2 = Minimap2()
        minimap2._print_stats(bam_file)

        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "samtools" in call_args
        assert "flagstat" in call_args

    @patch('circleseeker.external.minimap2.shutil.which')
    @patch('circleseeker.external.minimap2.Minimap2.align')
    def test_map_reads_legacy(self, mock_align, mock_which):
        """Test legacy map_reads method."""
        mock_which.side_effect = (
            lambda tool: f"/usr/bin/{tool}" if tool in {"minimap2", "samtools"} else None
        )

        reference = Path("reference.fasta")
        reads = Path("reads.fasta")
        output_sam = Path("output.sam")

        minimap2 = Minimap2()
        minimap2.map_reads(reference, reads, output_sam)

        mock_align.assert_called_once_with(reference, reads, output_sam, force_format="sam")


@pytest.fixture
def sample_reference_fasta(tmp_path):
    """Create a sample reference FASTA file."""
    fasta_file = tmp_path / "reference.fasta"
    fasta_content = """>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"""
    fasta_file.write_text(fasta_content)
    return fasta_file


@pytest.fixture
def sample_reads_fasta(tmp_path):
    """Create a sample reads FASTA file."""
    fasta_file = tmp_path / "reads.fasta"
    fasta_content = """>read1
ATCGATCGATCGATCG
>read2
GCTAGCTAGCTAGCTA
>read3
TTTTAAAAGGGGCCCC"""
    fasta_file.write_text(fasta_content)
    return fasta_file
