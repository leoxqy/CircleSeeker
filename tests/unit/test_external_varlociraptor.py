"""Tests for external Varlociraptor tools."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock, mock_open
import subprocess

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.varlociraptor import Varlociraptor


class TestVarlociraptor:
    """Test cases for Varlociraptor."""

    def test_varlociraptor_initialization(self):
        """Test Varlociraptor initialization."""
        varlociraptor = Varlociraptor()
        assert varlociraptor.tool_name == "varlociraptor"
        assert varlociraptor.threads == 1  # Default from base class

    def test_varlociraptor_initialization_with_threads(self):
        """Test Varlociraptor initialization with threads."""
        varlociraptor = Varlociraptor(threads=4)
        assert varlociraptor.threads == 4

    @patch.object(Varlociraptor, 'run')
    def test_estimate_alignment_properties(self, mock_run, tmp_path):
        """Test estimate_alignment_properties method."""
        # Mock return JSON output
        json_output = '{"insert_size": {"mean": 300, "std": 50}}'
        mock_run.return_value = (json_output, "")

        reference = tmp_path / "reference.fasta"
        bam = tmp_path / "input.bam"
        output_json = tmp_path / "output" / "alignprops.json"

        # Create input files
        reference.write_text(">chr1\nATCG\n")
        bam.touch()

        varlociraptor = Varlociraptor()
        varlociraptor.estimate_alignment_properties(reference, bam, output_json)

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "varlociraptor" in call_args
        assert "estimate" in call_args
        assert "alignment-properties" in call_args
        assert str(reference) in call_args
        assert "--bams" in call_args
        assert str(bam) in call_args

        # Verify capture_output parameter
        assert mock_run.call_args[1]['capture_output'] is True

        # Verify output directory is created and JSON is written
        assert output_json.parent.exists()
        assert output_json.exists()
        assert output_json.read_text() == json_output

    @patch.object(Varlociraptor, 'run')
    def test_preprocess_variants(self, mock_run, tmp_path):
        """Test preprocess_variants method."""
        reference = tmp_path / "reference.fasta"
        candidates_bcf = tmp_path / "candidates.bcf"
        alignprops_json = tmp_path / "alignprops.json"
        bam = tmp_path / "input.bam"
        output_obs_bcf = tmp_path / "output" / "observations.bcf"

        # Create input files
        reference.write_text(">chr1\nATCG\n")
        candidates_bcf.touch()
        alignprops_json.write_text('{"insert_size": {"mean": 300}}')
        bam.touch()

        varlociraptor = Varlociraptor()
        varlociraptor.preprocess_variants(
            reference=reference,
            candidates_bcf_sorted=candidates_bcf,
            alignprops_json=alignprops_json,
            bam=bam,
            output_obs_bcf=output_obs_bcf,
            max_depth=150
        )

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "varlociraptor" in call_args
        assert "preprocess" in call_args
        assert "variants" in call_args
        assert str(reference) in call_args
        assert "--candidates" in call_args
        assert str(candidates_bcf) in call_args
        assert "--alignment-properties" in call_args
        assert str(alignprops_json) in call_args
        assert "--max-depth" in call_args
        assert "150" in call_args
        assert "--bam" in call_args
        assert str(bam) in call_args
        assert "--output" in call_args
        assert str(output_obs_bcf) in call_args

        # Verify capture_output parameter
        assert mock_run.call_args[1]['capture_output'] is False

        # Verify output directory is created
        assert output_obs_bcf.parent.exists()

    @patch.object(Varlociraptor, 'run')
    def test_preprocess_variants_default_depth(self, mock_run, tmp_path):
        """Test preprocess_variants with default max_depth."""
        reference = tmp_path / "reference.fasta"
        candidates_bcf = tmp_path / "candidates.bcf"
        alignprops_json = tmp_path / "alignprops.json"
        bam = tmp_path / "input.bam"
        output_obs_bcf = tmp_path / "observations.bcf"

        reference.write_text(">chr1\nATCG\n")
        candidates_bcf.touch()
        alignprops_json.write_text('{}')
        bam.touch()

        varlociraptor = Varlociraptor()
        varlociraptor.preprocess_variants(
            reference=reference,
            candidates_bcf_sorted=candidates_bcf,
            alignprops_json=alignprops_json,
            bam=bam,
            output_obs_bcf=output_obs_bcf
        )

        call_args = mock_run.call_args[0][0]
        # Check default max_depth value
        assert "--max-depth" in call_args and "200" in call_args

    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_call_variants_generic(self, mock_file, mock_subprocess, tmp_path):
        """Test call_variants_generic method."""
        obs_bcf = tmp_path / "observations.bcf"
        scenario_yaml = tmp_path / "scenario.yaml"
        output_calls_bcf = tmp_path / "output" / "calls.bcf"

        # Create input files
        obs_bcf.touch()
        scenario_yaml.write_text("events: [PRESENT]")

        mock_subprocess.return_value = MagicMock(returncode=0)

        varlociraptor = Varlociraptor()
        varlociraptor.call_variants_generic(
            obs_bcf_sorted=obs_bcf,
            sample_name="test_sample",
            scenario_yaml=scenario_yaml,
            output_calls_bcf=output_calls_bcf
        )

        # Verify subprocess was called with correct command
        mock_subprocess.assert_called_once()
        call_args = mock_subprocess.call_args[0][0]
        assert "varlociraptor" in call_args
        assert "call" in call_args
        assert "variants" in call_args
        assert "generic" in call_args
        assert "--obs" in call_args
        assert f"test_sample={obs_bcf}" in call_args
        assert "--scenario" in call_args
        assert str(scenario_yaml) in call_args

        # Verify file was opened for writing
        mock_file.assert_called_once_with(output_calls_bcf, "wb")

        # Verify output directory is created
        assert output_calls_bcf.parent.exists()

    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_call_variants_generic_failure(self, mock_file, mock_subprocess, tmp_path):
        """Test call_variants_generic with subprocess failure."""
        from subprocess import CalledProcessError

        obs_bcf = tmp_path / "observations.bcf"
        scenario_yaml = tmp_path / "scenario.yaml"
        output_calls_bcf = tmp_path / "calls.bcf"

        obs_bcf.touch()
        scenario_yaml.write_text("events: [PRESENT]")

        # Mock subprocess failure
        mock_subprocess.side_effect = CalledProcessError(1, "varlociraptor", stderr=b"Error message")

        varlociraptor = Varlociraptor()

        with pytest.raises(Exception):  # Should raise PipelineError
            varlociraptor.call_variants_generic(
                obs_bcf_sorted=obs_bcf,
                sample_name="test_sample",
                scenario_yaml=scenario_yaml,
                output_calls_bcf=output_calls_bcf
            )

    @patch('subprocess.Popen')
    def test_filter_calls_fdr_local_smart(self, mock_popen, tmp_path):
        """Test filter_calls_fdr_local_smart method."""
        input_calls_bcf = tmp_path / "calls.bcf"
        output_calls_fdr_bcf = tmp_path / "output" / "calls_filtered.bcf"

        input_calls_bcf.touch()

        # Mock successful pipeline processes
        mock_p1 = MagicMock()
        mock_p1.wait.return_value = 0
        mock_p1.stdout = MagicMock()

        mock_p2 = MagicMock()
        mock_p2.wait.return_value = 0
        mock_p2.stdout = MagicMock()

        mock_p3 = MagicMock()
        mock_p3.wait.return_value = 0

        mock_popen.side_effect = [mock_p1, mock_p2, mock_p3]

        varlociraptor = Varlociraptor()
        varlociraptor.filter_calls_fdr_local_smart(
            input_calls_bcf=input_calls_bcf,
            output_calls_fdr_bcf=output_calls_fdr_bcf,
            fdr=1.0,
            memory_limit="8G"
        )

        # Verify all three processes were created
        assert mock_popen.call_count == 3

        # Check first command (filter-calls)
        first_call = mock_popen.call_args_list[0][0][0]
        assert "varlociraptor" in first_call
        assert "filter-calls" in first_call
        assert "control-fdr" in first_call
        assert "--fdr" in first_call
        assert "1.0" in first_call
        assert str(input_calls_bcf) in first_call

        # Check second command (decode-phred)
        second_call = mock_popen.call_args_list[1][0][0]
        assert "varlociraptor" in second_call
        assert "decode-phred" in second_call

        # Check third command (bcftools sort)
        third_call = mock_popen.call_args_list[2][0][0]
        assert "bcftools" in third_call
        assert "sort" in third_call
        assert "-m" in third_call
        assert "8G" in third_call
        assert str(output_calls_fdr_bcf) in third_call

        # Verify output directory is created
        assert output_calls_fdr_bcf.parent.exists()

    @patch('subprocess.Popen')
    def test_filter_calls_fdr_defaults(self, mock_popen, tmp_path):
        """Test filter_calls_fdr_local_smart with default parameters."""
        input_calls_bcf = tmp_path / "calls.bcf"
        output_calls_fdr_bcf = tmp_path / "calls_filtered.bcf"

        input_calls_bcf.touch()

        # Mock successful processes
        mock_p1 = MagicMock()
        mock_p1.wait.return_value = 0
        mock_p1.stdout = MagicMock()

        mock_p2 = MagicMock()
        mock_p2.wait.return_value = 0
        mock_p2.stdout = MagicMock()

        mock_p3 = MagicMock()
        mock_p3.wait.return_value = 0

        mock_popen.side_effect = [mock_p1, mock_p2, mock_p3]

        varlociraptor = Varlociraptor()
        varlociraptor.filter_calls_fdr_local_smart(
            input_calls_bcf=input_calls_bcf,
            output_calls_fdr_bcf=output_calls_fdr_bcf
        )

        # Check default values in first command
        first_call = mock_popen.call_args_list[0][0][0]
        assert "--fdr" in first_call and "0.2" in first_call

        # Check default memory in third command
        third_call = mock_popen.call_args_list[2][0][0]
        assert "-m" in third_call and "4G" in third_call

    @patch('subprocess.Popen')
    def test_filter_calls_fdr_pipeline_failure(self, mock_popen, tmp_path):
        """Test filter_calls_fdr_local_smart with pipeline failure."""
        input_calls_bcf = tmp_path / "calls.bcf"
        output_calls_fdr_bcf = tmp_path / "calls_filtered.bcf"

        input_calls_bcf.touch()

        # Mock failed first process
        mock_p1 = MagicMock()
        mock_p1.wait.return_value = 1
        mock_p1.communicate.return_value = (None, b"Filter failed")
        mock_p1.stdout = MagicMock()

        mock_p2 = MagicMock()
        mock_p2.wait.return_value = 0
        mock_p2.stdout = MagicMock()

        mock_p3 = MagicMock()
        mock_p3.wait.return_value = 0

        mock_popen.side_effect = [mock_p1, mock_p2, mock_p3]

        varlociraptor = Varlociraptor()

        with pytest.raises(Exception):  # Should raise PipelineError
            varlociraptor.filter_calls_fdr_local_smart(
                input_calls_bcf=input_calls_bcf,
                output_calls_fdr_bcf=output_calls_fdr_bcf
            )


@pytest.fixture
def sample_bcf_file(tmp_path):
    """Create a sample BCF file for testing."""
    bcf_file = tmp_path / "sample.bcf"
    # Create a dummy BCF file (in real scenario this would be binary)
    bcf_file.write_bytes(b"BCF\x02\x02")  # Minimal BCF header
    return bcf_file


@pytest.fixture
def sample_scenario_yaml(tmp_path):
    """Create a sample scenario YAML file for testing."""
    yaml_file = tmp_path / "scenario.yaml"
    yaml_content = """events:
  PRESENT:
    prior: 0.1
    description: "Present variant"
"""
    yaml_file.write_text(yaml_content)
    return yaml_file


@pytest.fixture
def sample_alignprops_json(tmp_path):
    """Create a sample alignment properties JSON file for testing."""
    json_file = tmp_path / "alignprops.json"
    json_content = """{
  "insert_size": {
    "mean": 300.5,
    "std": 50.2
  },
  "max_del_cigar_len": 1000,
  "max_ins_cigar_len": 500
}"""
    json_file.write_text(json_content)
    return json_file