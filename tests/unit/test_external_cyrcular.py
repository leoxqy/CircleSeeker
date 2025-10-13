"""Tests for external Cyrcular tools."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.cyrcular import Cyrcular


class TestCyrcular:
    """Test cases for Cyrcular."""

    def test_cyrcular_initialization(self):
        """Test Cyrcular initialization."""
        cyrcular = Cyrcular()
        assert cyrcular.tool_name == "cyrcular"
        assert cyrcular.threads == 1  # Default from base class

    def test_cyrcular_initialization_with_threads(self):
        """Test Cyrcular initialization with threads."""
        cyrcular = Cyrcular(threads=4)
        assert cyrcular.threads == 4

    @patch.object(Cyrcular, 'run')
    def test_graph_breakends(self, mock_run, tmp_path):
        """Test graph_breakends method."""
        bam = tmp_path / "input.bam"
        reference = tmp_path / "reference.fasta"
        output_candidates = tmp_path / "output" / "candidates.bcf"
        output_graph = tmp_path / "output" / "graph.json"
        dot_dir = tmp_path / "dot_output"

        # Create input files
        bam.touch()
        reference.write_text(">chr1\nATCG\n")

        cyrcular = Cyrcular()
        cyrcular.graph_breakends(
            bam=bam,
            reference=reference,
            output_candidates=output_candidates,
            output_graph=output_graph,
            dot_dir=dot_dir,
            min_read_depth=5,
            min_split_reads=4,
            max_paths_per_component=15,
            max_deletion_length=500,
            threads=2
        )

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "cyrcular" in call_args
        assert "graph" in call_args
        assert "breakends" in call_args
        assert str(bam) in call_args
        assert "--reference" in call_args
        assert str(reference) in call_args
        assert "--min-read-depth" in call_args
        assert "5" in call_args
        assert "--min-split-reads" in call_args
        assert "4" in call_args
        assert "--max-paths-per-component" in call_args
        assert "15" in call_args
        assert "--max-deletion-length" in call_args
        assert "500" in call_args
        assert "-t" in call_args
        assert "2" in call_args
        assert "--output" in call_args
        assert str(output_candidates) in call_args
        assert "--graph" in call_args
        assert str(output_graph) in call_args
        assert "--dot" in call_args

        # Verify capture_output parameter
        assert mock_run.call_args[1]['capture_output'] is False

        # Verify directories are created
        assert output_candidates.parent.exists()
        assert dot_dir.exists()

    @patch.object(Cyrcular, 'run')
    def test_graph_breakends_defaults(self, mock_run, tmp_path):
        """Test graph_breakends with default parameters."""
        bam = tmp_path / "input.bam"
        reference = tmp_path / "reference.fasta"
        output_candidates = tmp_path / "candidates.bcf"
        output_graph = tmp_path / "graph.json"
        dot_dir = tmp_path / "dot"

        bam.touch()
        reference.write_text(">chr1\nATCG\n")

        cyrcular = Cyrcular()
        cyrcular.graph_breakends(
            bam=bam,
            reference=reference,
            output_candidates=output_candidates,
            output_graph=output_graph,
            dot_dir=dot_dir
        )

        call_args = mock_run.call_args[0][0]
        # Check default values
        assert "--min-read-depth" in call_args and "2" in call_args
        assert "--min-split-reads" in call_args and "2" in call_args
        assert "--max-paths-per-component" in call_args and "20" in call_args
        assert "--max-deletion-length" in call_args and "1000" in call_args
        assert "-t" in call_args and "1" in call_args

    @patch.object(Cyrcular, 'run')
    def test_graph_annotate(self, mock_run, tmp_path):
        """Test graph_annotate method."""
        reference = tmp_path / "reference.fasta"
        gene_annotation = tmp_path / "genes.gff3.gz"
        regulatory_annotation = tmp_path / "regulatory.gff3.gz"
        graph_input = tmp_path / "input_graph.json"
        output_graph = tmp_path / "output" / "annotated_graph.json"

        # Create input files
        reference.write_text(">chr1\nATCG\n")
        gene_annotation.touch()
        regulatory_annotation.touch()
        graph_input.write_text("{}")

        cyrcular = Cyrcular()
        cyrcular.graph_annotate(
            reference=reference,
            gene_annotation_gff_gz=gene_annotation,
            regulatory_annotation_gff_gz=regulatory_annotation,
            graph_input=graph_input,
            output_graph=output_graph
        )

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "cyrcular" in call_args
        assert "graph" in call_args
        assert "annotate" in call_args
        assert "--reference" in call_args
        assert str(reference) in call_args
        assert "--gene-annotation" in call_args
        assert str(gene_annotation) in call_args
        assert "--regulatory-annotation" in call_args
        assert str(regulatory_annotation) in call_args
        assert "--output" in call_args
        assert str(output_graph) in call_args
        assert str(graph_input) in call_args

        # Verify capture_output parameter
        assert mock_run.call_args[1]['capture_output'] is False

        # Verify output directory is created
        assert output_graph.parent.exists()

    @patch.object(Cyrcular, 'run')
    def test_graph_table(self, mock_run, tmp_path):
        """Test graph_table method."""
        annotated_graph = tmp_path / "annotated_graph.json"
        calls_bcf = tmp_path / "calls.bcf"
        reference = tmp_path / "reference.fasta"
        circle_table = tmp_path / "output" / "circles.tsv"
        segment_tables_dir = tmp_path / "segments"

        # Create input files
        annotated_graph.write_text("{}")
        calls_bcf.touch()
        reference.write_text(">chr1\nATCG\n")

        cyrcular = Cyrcular()
        cyrcular.graph_table(
            annotated_graph=annotated_graph,
            calls_bcf=calls_bcf,
            reference=reference,
            circle_table=circle_table,
            segment_tables_dir=segment_tables_dir
        )

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "cyrcular" in call_args
        assert "graph" in call_args
        assert "table" in call_args
        assert str(annotated_graph) in call_args
        assert str(calls_bcf) in call_args
        assert "--reference" in call_args
        assert str(reference) in call_args
        assert "--circle-table" in call_args
        assert str(circle_table) in call_args
        assert "--segment-tables" in call_args

        # Verify capture_output parameter
        assert mock_run.call_args[1]['capture_output'] is False

        # Verify directories are created
        assert circle_table.parent.exists()
        assert segment_tables_dir.exists()

    @patch.object(Cyrcular, 'run')
    def test_complete_workflow(self, mock_run, tmp_path):
        """Test complete cyrcular workflow."""
        # Setup files
        bam = tmp_path / "input.bam"
        reference = tmp_path / "reference.fasta"
        gene_annotation = tmp_path / "genes.gff3.gz"
        regulatory_annotation = tmp_path / "regulatory.gff3.gz"

        # Intermediate files
        candidates = tmp_path / "candidates.bcf"
        graph = tmp_path / "graph.json"
        annotated_graph = tmp_path / "annotated_graph.json"
        dot_dir = tmp_path / "dot"

        # Output files
        circle_table = tmp_path / "circles.tsv"
        segment_tables_dir = tmp_path / "segments"

        # Create input files
        bam.touch()
        reference.write_text(">chr1\nATCG\n")
        gene_annotation.touch()
        regulatory_annotation.touch()

        cyrcular = Cyrcular()

        # Step 1: Breakend detection
        cyrcular.graph_breakends(
            bam=bam,
            reference=reference,
            output_candidates=candidates,
            output_graph=graph,
            dot_dir=dot_dir
        )

        # Step 2: Annotation
        cyrcular.graph_annotate(
            reference=reference,
            gene_annotation_gff_gz=gene_annotation,
            regulatory_annotation_gff_gz=regulatory_annotation,
            graph_input=graph,
            output_graph=annotated_graph
        )

        # Step 3: Table generation
        cyrcular.graph_table(
            annotated_graph=annotated_graph,
            calls_bcf=candidates,
            reference=reference,
            circle_table=circle_table,
            segment_tables_dir=segment_tables_dir
        )

        # Verify all three operations were called
        assert mock_run.call_count == 3

        # Check first call (breakends)
        first_call = mock_run.call_args_list[0][0][0]
        assert "breakends" in first_call

        # Check second call (annotate)
        second_call = mock_run.call_args_list[1][0][0]
        assert "annotate" in second_call

        # Check third call (table)
        third_call = mock_run.call_args_list[2][0][0]
        assert "table" in third_call

    @patch.object(Cyrcular, 'run')
    def test_error_handling(self, mock_run, tmp_path):
        """Test error handling when cyrcular operations fail."""
        mock_run.side_effect = Exception("Cyrcular command failed")

        bam = tmp_path / "input.bam"
        reference = tmp_path / "reference.fasta"
        output_candidates = tmp_path / "candidates.bcf"
        output_graph = tmp_path / "graph.json"
        dot_dir = tmp_path / "dot"

        bam.touch()
        reference.write_text(">chr1\nATCG\n")

        cyrcular = Cyrcular()

        with pytest.raises(Exception):
            cyrcular.graph_breakends(
                bam=bam,
                reference=reference,
                output_candidates=output_candidates,
                output_graph=output_graph,
                dot_dir=dot_dir
            )


@pytest.fixture
def sample_bam_file(tmp_path):
    """Create a sample BAM file for testing."""
    bam_file = tmp_path / "sample.bam"
    # Create a dummy BAM file (in real scenario this would be binary)
    bam_file.write_bytes(b"BAM\x01")  # Minimal BAM header
    return bam_file


@pytest.fixture
def sample_reference_fasta(tmp_path):
    """Create a sample reference FASTA file for testing."""
    ref_file = tmp_path / "reference.fasta"
    ref_content = """>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"""
    ref_file.write_text(ref_content)
    return ref_file


@pytest.fixture
def sample_gff_file(tmp_path):
    """Create a sample GFF3 file for testing."""
    gff_file = tmp_path / "annotations.gff3.gz"
    # Create a dummy compressed GFF file
    gff_file.write_bytes(b"\x1f\x8b")  # GZip magic bytes
    return gff_file