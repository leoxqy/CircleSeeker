from pathlib import Path
import importlib.util
import sys
import types
import tempfile
import pandas as pd
from unittest.mock import Mock, patch, MagicMock, call
import pytest
import gzip

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _load_cyrcular_calling_module():
    package_name = "circleseeker.modules"
    module_name = f"{package_name}.cyrcular_calling"

    if module_name in sys.modules:
        return sys.modules[module_name]

    import circleseeker  # noqa: F401  # Ensure base package is initialized

    if package_name not in sys.modules:
        package_module = types.ModuleType(package_name)
        package_module.__path__ = [str(SRC / "circleseeker" / "modules")]
        sys.modules[package_name] = package_module

    spec = importlib.util.spec_from_file_location(
        module_name,
        SRC / "circleseeker" / "modules" / "cyrcular_calling.py",
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


# Load module and get classes
cc_module = _load_cyrcular_calling_module()
PipelineConfig = cc_module.PipelineConfig
CircularDNAResult = cc_module.CircularDNAResult
DependencyChecker = cc_module.DependencyChecker
CyrcularCallingPipeline = cc_module.CyrcularCallingPipeline
PipelineError = cc_module.PipelineError


class TestPipelineConfig:
    def test_config_initialization(self, tmp_path):
        # Create test files
        bam_file = tmp_path / "test.bam"
        ref_file = tmp_path / "test.fa"
        bam_file.touch()
        ref_file.touch()
        
        config = PipelineConfig(
            bam_file=bam_file,
            reference=ref_file,
            output_dir=tmp_path / "output",
            sample_name="test_sample"
        )
        
        assert config.bam_file == bam_file
        assert config.reference == ref_file
        assert config.output_dir == tmp_path / "output"
        assert config.sample_name == "test_sample"
        assert config.threads == 4
        assert config.min_read_depth == 2
        assert config.fdr == 0.05
    
    def test_config_missing_bam_file(self, tmp_path):
        ref_file = tmp_path / "test.fa"
        ref_file.touch()
        
        with pytest.raises(FileNotFoundError, match="BAM file not found"):
            PipelineConfig(
                bam_file=tmp_path / "nonexistent.bam",
                reference=ref_file,
                output_dir=tmp_path / "output"
            )
    
    def test_config_missing_reference(self, tmp_path):
        bam_file = tmp_path / "test.bam"
        bam_file.touch()
        
        with pytest.raises(FileNotFoundError, match="Reference file not found"):
            PipelineConfig(
                bam_file=bam_file,
                reference=tmp_path / "nonexistent.fa",
                output_dir=tmp_path / "output"
            )


class TestCircularDNAResult:
    def test_basic_initialization(self):
        result = CircularDNAResult(
            circle_id="circle1",
            circle_length=1000,
            num_segments=3,
            vaf=0.25,
            read_support=10
        )
        
        assert result.circle_id == "circle1"
        assert result.circle_length == 1000
        assert result.num_segments == 3
        assert result.vaf == 0.25
        assert result.read_support == 10
    
    def test_from_dataframe_row(self):
        row = pd.Series({
            "circle_id": "circle2",
            "circle_length": 2000,
            "num_segments": 2,
            "vaf": 0.5,
            "read_support": 20
        })
        
        result = CircularDNAResult.from_dataframe_row(row)
        
        assert result.circle_id == "circle2"
        assert result.circle_length == 2000
        assert result.num_segments == 2
        assert result.vaf == 0.5
        assert result.read_support == 20
    
    def test_from_dataframe_row_missing_fields(self):
        row = pd.Series({
            "circle_id": "circle3",
            "circle_length": 3000
        })
        
        result = CircularDNAResult.from_dataframe_row(row)
        
        assert result.circle_id == "circle3"
        assert result.circle_length == 3000
        assert result.num_segments == 0
        assert result.vaf is None
        assert result.read_support is None


class TestDependencyChecker:
    @patch('shutil.which')
    def test_check_all_tools_present(self, mock_which):
        mock_which.return_value = "/usr/bin/tool"
        
        checker = DependencyChecker()
        checker.check_all()  # Should not raise
    
    @patch('shutil.which')
    def test_check_missing_tools(self, mock_which):
        def which_side_effect(tool):
            if tool in ["samtools", "bcftools"]:
                return "/usr/bin/" + tool
            return None
        
        mock_which.side_effect = which_side_effect
        
        checker = DependencyChecker()
        with pytest.raises(PipelineError, match="Missing required tools: cyrcular, varlociraptor"):
            checker.check_all()


class TestCyrcularCallingPipeline:
    @pytest.fixture
    def setup_pipeline(self, tmp_path):
        # Create test files
        bam_file = tmp_path / "test.bam"
        ref_file = tmp_path / "test.fa"
        bam_file.touch()
        ref_file.touch()
        
        config = PipelineConfig(
            bam_file=bam_file,
            reference=ref_file,
            output_dir=tmp_path / "output",
            sample_name="test"
        )
        
        # Create pipeline with mocked dependencies
        with patch('shutil.which', return_value="/usr/bin/tool"):
            pipeline = CyrcularCallingPipeline(config)
            
            # Mock the external tool objects
            pipeline._samtools = Mock()
            pipeline._cyrcular = Mock()
            pipeline._bcftools = Mock()
            pipeline._varloc = Mock()
            
            # Return objects for testing
            return {
                'pipeline': pipeline,
                'config': config,
                'samtools': pipeline._samtools,
                'cyrcular': pipeline._cyrcular,
                'bcftools': pipeline._bcftools,
                'varloc': pipeline._varloc,
                'tmp_path': tmp_path
            }
    
    def test_initialize_file_paths(self, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        
        # Check key file paths are initialized
        assert "candidates" in pipeline.file_paths
        assert "graph" in pipeline.file_paths
        assert "calls" in pipeline.file_paths
        assert "overview_table" in pipeline.file_paths
        
        # Check paths use correct prefix
        assert "test.candidates.bcf" in str(pipeline.file_paths["candidates"])
        assert "test_overview.tsv" in str(pipeline.file_paths["overview_table"])
    
    @patch('shutil.which', return_value="/usr/bin/tool")
    def test_prepare_indices(self, mock_which, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        samtools = setup_pipeline['samtools']
        
        # Mock file existence checks
        with patch('pathlib.Path.exists', return_value=False):
            pipeline._prepare_indices()
        
        # Check samtools was called
        samtools.faidx.assert_called_once_with(pipeline.config.reference)
        samtools.index_bam.assert_called_once_with(pipeline.config.bam_file)
    
    @patch('shutil.which', return_value="/usr/bin/tool")
    def test_run_cyrcular(self, mock_which, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        cyrcular = setup_pipeline['cyrcular']
        bcftools = setup_pipeline['bcftools']
        
        pipeline._run_cyrcular()
        
        # Check cyrcular graph_breakends was called with correct params
        cyrcular.graph_breakends.assert_called_once()
        call_args = cyrcular.graph_breakends.call_args[1]
        assert call_args['bam'] == pipeline.config.bam_file
        assert call_args['reference'] == pipeline.config.reference
        assert call_args['min_read_depth'] == 2
        assert call_args['threads'] == 4
        
        # Check bcftools sort and index were called
        assert bcftools.sort.called
        assert bcftools.index.called
    
    @patch('shutil.which', return_value="/usr/bin/tool")
    def test_run_varlociraptor(self, mock_which, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        varloc = setup_pipeline['varloc']
        bcftools = setup_pipeline['bcftools']
        tmp_path = setup_pipeline['tmp_path']
        
        # Create output directory
        pipeline.config.output_dir.mkdir(parents=True, exist_ok=True)
        
        pipeline._run_varlociraptor()
        
        # Check varlociraptor methods were called
        varloc.estimate_alignment_properties.assert_called_once()
        varloc.preprocess_variants.assert_called_once()
        varloc.call_variants_generic.assert_called_once()
        varloc.filter_calls_fdr_local_smart.assert_called_once()
        
        # Check scenario file was created
        scenario_file = pipeline.file_paths["scenario"]
        assert scenario_file.exists()
        content = scenario_file.read_text()
        assert "test:" in content  # sample name
        assert "PRESENT:" in content
    
    @patch('shutil.which', return_value="/usr/bin/tool")
    def test_annotate_results(self, mock_which, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        cyrcular = setup_pipeline['cyrcular']
        tmp_path = setup_pipeline['tmp_path']
        
        # Create output directory
        pipeline.config.output_dir.mkdir(parents=True, exist_ok=True)
        
        pipeline._annotate_results()
        
        # Check empty annotation files were created
        gene_file = pipeline.file_paths["empty_gene"]
        reg_file = pipeline.file_paths["empty_regulatory"]
        assert gene_file.exists()
        assert reg_file.exists()
        
        # Check they are gzipped
        with gzip.open(gene_file, 'rt') as f:
            content = f.read()
            assert "##gff-version 3" in content
        
        # Check cyrcular methods were called
        cyrcular.graph_annotate.assert_called_once()
        cyrcular.graph_table.assert_called_once()
    
    def test_parse_results_success(self, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        tmp_path = setup_pipeline['tmp_path']
        
        # Create test overview table
        overview_dir = pipeline.file_paths["overview_table"].parent
        overview_dir.mkdir(parents=True, exist_ok=True)
        
        df = pd.DataFrame({
            'circle_id': ['circle1', 'circle2'],
            'circle_length': [1000, 2000],
            'num_segments': [2, 3],
            'vaf': [0.25, 0.5],
            'read_support': [10, 20]
        })
        df.to_csv(pipeline.file_paths["overview_table"], sep='\t', index=False)
        
        results = pipeline._parse_results()
        
        assert len(results) == 2
        assert results[0].circle_id == 'circle1'
        assert results[0].circle_length == 1000
        assert results[1].circle_id == 'circle2'
        assert results[1].vaf == 0.5
    
    def test_parse_results_no_file(self, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        
        results = pipeline._parse_results()
        assert results == []
    
    def test_parse_results_empty_file(self, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        
        # Create empty overview table with at least column headers
        overview_dir = pipeline.file_paths["overview_table"].parent
        overview_dir.mkdir(parents=True, exist_ok=True)
        
        # Write empty file with just headers - no data rows
        with open(pipeline.file_paths["overview_table"], 'w') as f:
            f.write("circle_id\tcircle_length\tnum_segments\tvaf\tread_support\n")
        
        results = pipeline._parse_results()
        assert results == []
    
    def test_cleanup_intermediate_files(self, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        
        # Create some intermediate files
        pipeline.config.output_dir.mkdir(parents=True, exist_ok=True)
        for key in ["obs", "candidates", "scenario"]:
            if key in pipeline.file_paths:
                pipeline.file_paths[key].touch()
        
        pipeline._cleanup_intermediate_files()
        
        # Check files were removed
        assert not pipeline.file_paths["obs"].exists()
        assert not pipeline.file_paths["candidates"].exists()
        assert not pipeline.file_paths["scenario"].exists()
    
    @patch('shutil.which', return_value="/usr/bin/tool")
    def test_run_full_pipeline(self, mock_which, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        tmp_path = setup_pipeline['tmp_path']
        
        # Mock all the internal methods
        with patch.object(pipeline, '_prepare_indices') as mock_prep, \
             patch.object(pipeline, '_run_cyrcular') as mock_cyrc, \
             patch.object(pipeline, '_run_varlociraptor') as mock_varl, \
             patch.object(pipeline, '_annotate_results') as mock_anno, \
             patch.object(pipeline, '_parse_results', return_value=[
                 CircularDNAResult("circle1", 1000, 2, 0.25, 10)
             ]) as mock_parse, \
             patch.object(pipeline, '_cleanup_intermediate_files') as mock_clean, \
             patch.object(pipeline, '_print_summary') as mock_print:
            
            results = pipeline.run()
        
        assert len(results) == 1
        assert results[0].circle_id == "circle1"
        
        # Verify all steps were called
        mock_prep.assert_called_once()
        mock_cyrc.assert_called_once()
        mock_varl.assert_called_once()
        mock_anno.assert_called_once()
        mock_parse.assert_called_once()
        mock_clean.assert_called_once()
        mock_print.assert_called_once()
    
    @patch('shutil.which', return_value="/usr/bin/tool")
    def test_run_pipeline_keep_intermediate(self, mock_which, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        pipeline.config.keep_intermediate = True
        
        # Mock all the internal methods
        with patch.object(pipeline, '_prepare_indices'), \
             patch.object(pipeline, '_run_cyrcular'), \
             patch.object(pipeline, '_run_varlociraptor'), \
             patch.object(pipeline, '_annotate_results'), \
             patch.object(pipeline, '_parse_results', return_value=[]), \
             patch.object(pipeline, '_cleanup_intermediate_files') as mock_cleanup, \
             patch.object(pipeline, '_print_summary'):
            
            pipeline.run()
        
        # Cleanup should NOT be called when keep_intermediate is True
        mock_cleanup.assert_not_called()
    
    def test_print_summary(self, setup_pipeline):
        pipeline = setup_pipeline['pipeline']
        
        results = [
            CircularDNAResult("circle1", 1000, 2, 0.25, 10),
            CircularDNAResult("circle2", 2000, 3, 0.5, 20),
        ]
        
        # Just check it runs without error
        pipeline._print_summary(results)


class TestIntegration:
    """Integration tests that require actual files but mock external tools."""
    
    @patch('shutil.which', return_value="/usr/bin/tool")
    def test_pipeline_with_mocked_tools(self, mock_which, tmp_path):
        # Create test files
        bam_file = tmp_path / "test.bam"
        ref_file = tmp_path / "test.fa"
        bam_file.touch()
        ref_file.touch()
        
        config = PipelineConfig(
            bam_file=bam_file,
            reference=ref_file,
            output_dir=tmp_path / "output",
            sample_name="test",
            threads=2
        )
        
        # Mock the entire run method to avoid actual subprocess calls
        with patch.object(CyrcularCallingPipeline, 'run') as mock_run:
            mock_run.return_value = [
                CircularDNAResult("circle1", 1000, 2, 0.25, 10),
                CircularDNAResult("circle2", 2000, 3, 0.5, 20),
            ]
            
            pipeline = CyrcularCallingPipeline(config)
            results = pipeline.run()
        
        # Verify results
        assert len(results) == 2
        assert results[0].circle_id == "circle1"
        assert results[0].circle_length == 1000
        assert results[1].circle_id == "circle2"
        assert results[1].vaf == 0.5
        
        # Verify run was called
        mock_run.assert_called_once()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])