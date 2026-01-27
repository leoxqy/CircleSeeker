from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import pytest

from circleseeker.core.steps.contracts import (
    ArtifactSpec,
    StepContract,
    STEP_CONTRACTS,
    get_step_contract,
)


# ---------------------------------------------------------------------------
# TestArtifactSpec
# ---------------------------------------------------------------------------
class TestArtifactSpec:
    """Tests for the ArtifactSpec dataclass."""

    def test_create_with_required_fields(self):
        spec = ArtifactSpec(name="test", base="output", template="test.csv")
        assert spec.name == "test"
        assert spec.base == "output"
        assert spec.template == "test.csv"

    def test_defaults_kind_file_required_false(self):
        spec = ArtifactSpec(name="x", base="output", template="x.txt")
        assert spec.kind == "file"
        assert spec.required is False

    def test_required_columns_tuple_stored(self):
        cols = ("chr", "start0", "end0")
        spec = ArtifactSpec(
            name="t",
            base="output",
            template="t.csv",
            required_columns=cols,
        )
        assert spec.required_columns == ("chr", "start0", "end0")
        assert isinstance(spec.required_columns, tuple)


# ---------------------------------------------------------------------------
# TestStepContract
# ---------------------------------------------------------------------------
class TestStepContract:
    """Tests for the StepContract dataclass."""

    def test_create_with_fields(self):
        inp = ArtifactSpec(name="in1", base="output", template="in1.csv")
        out = ArtifactSpec(name="out1", base="output", template="out1.csv")
        contract = StepContract(step_name="my_step", inputs=(inp,), outputs=(out,))
        assert contract.step_name == "my_step"
        assert len(contract.inputs) == 1
        assert len(contract.outputs) == 1

    def test_default_empty_tuples(self):
        contract = StepContract(step_name="empty")
        assert contract.inputs == ()
        assert contract.outputs == ()


# ---------------------------------------------------------------------------
# TestStepContracts
# ---------------------------------------------------------------------------
class TestStepContracts:
    """Tests for the STEP_CONTRACTS mapping."""

    def test_step_contracts_is_dict(self):
        assert isinstance(STEP_CONTRACTS, dict)

    def test_contains_expected_step_names(self):
        expected = [
            "tidehunter",
            "um_classify",
            "cecc_build",
            "ecc_dedup",
            "ecc_unify",
            "ecc_packager",
        ]
        for name in expected:
            assert name in STEP_CONTRACTS, f"Missing step: {name}"

    def test_each_value_is_step_contract(self):
        for key, value in STEP_CONTRACTS.items():
            assert isinstance(value, StepContract), (
                f"STEP_CONTRACTS[{key!r}] is {type(value).__name__}, expected StepContract"
            )


# ---------------------------------------------------------------------------
# TestGetStepContract
# ---------------------------------------------------------------------------
class TestGetStepContract:
    """Tests for get_step_contract()."""

    def test_known_step_returns_contract(self):
        contract = get_step_contract("tidehunter")
        assert isinstance(contract, StepContract)
        assert contract.step_name == "tidehunter"

    def test_unknown_step_returns_none(self):
        result = get_step_contract("nonexistent_step_xyz")
        assert result is None
