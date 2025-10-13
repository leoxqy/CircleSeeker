"""Tests for utils progress module."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock
from typing import List

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.utils.progress import iter_progress


class TestIterProgress:
    """Test cases for iter_progress function."""

    def test_iter_progress_basic_functionality(self):
        """Test iter_progress basic functionality without tqdm."""
        test_data = [1, 2, 3, 4, 5]

        # Should work even without tqdm
        result = list(iter_progress(test_data))
        assert result == test_data

    def test_iter_progress_with_enabled_false(self):
        """Test iter_progress with enabled=False."""
        test_data = ['a', 'b', 'c']

        result = list(iter_progress(test_data, enabled=False))
        assert result == test_data

    def test_iter_progress_returns_iterator(self):
        """Test that iter_progress returns an iterator."""
        test_data = [1, 2, 3]

        progress_iter = iter_progress(test_data)

        # Should be an iterator
        assert hasattr(progress_iter, '__iter__')
        assert hasattr(progress_iter, '__next__')

        # Should be consumable
        result = list(progress_iter)
        assert result == test_data

    def test_iter_progress_with_different_iterables(self):
        """Test iter_progress with different types of iterables."""
        # Test with list
        list_data = [1, 2, 3]
        assert list(iter_progress(list_data)) == list_data

        # Test with tuple
        tuple_data = (4, 5, 6)
        assert list(iter_progress(tuple_data)) == list(tuple_data)

        # Test with range
        range_data = range(3)
        assert list(iter_progress(range_data)) == list(range_data)

        # Test with string (iterable)
        string_data = "abc"
        assert list(iter_progress(string_data)) == list(string_data)

    def test_iter_progress_empty_iterable(self):
        """Test iter_progress with empty iterable."""
        empty_list = []
        result = list(iter_progress(empty_list))
        assert result == []

    def test_iter_progress_single_item(self):
        """Test iter_progress with single item."""
        single_item = [42]
        result = list(iter_progress(single_item))
        assert result == [42]

    @patch('circleseeker.utils.progress.tqdm')
    def test_iter_progress_with_tqdm_available(self, mock_tqdm):
        """Test iter_progress when tqdm is available."""
        test_data = [1, 2, 3, 4]
        mock_tqdm_instance = MagicMock()
        mock_tqdm_instance.__iter__ = MagicMock(return_value=iter(test_data))
        mock_tqdm.return_value = mock_tqdm_instance

        result = list(iter_progress(test_data, total=4, desc="Testing"))

        # Should call tqdm with correct parameters
        mock_tqdm.assert_called_once_with(test_data, total=4, desc="Testing")
        assert result == test_data

    @patch('builtins.__import__')
    def test_iter_progress_tqdm_import_error(self, mock_import):
        """Test iter_progress when tqdm import fails."""
        # Mock import to raise ImportError for tqdm
        def import_side_effect(name, *args, **kwargs):
            if name == 'tqdm':
                raise ImportError("No module named 'tqdm'")
            return __import__(name, *args, **kwargs)

        mock_import.side_effect = import_side_effect

        test_data = [1, 2, 3]
        result = list(iter_progress(test_data))

        # Should fall back to plain iterator
        assert result == test_data

    @patch('circleseeker.utils.progress.tqdm')
    def test_iter_progress_tqdm_runtime_error(self, mock_tqdm):
        """Test iter_progress when tqdm raises runtime error."""
        test_data = [1, 2, 3]
        mock_tqdm.side_effect = RuntimeError("tqdm error")

        result = list(iter_progress(test_data))

        # Should fall back to plain iterator
        assert result == test_data

    @patch('circleseeker.utils.progress.tqdm')
    def test_iter_progress_with_total_parameter(self, mock_tqdm):
        """Test iter_progress with total parameter."""
        test_data = [1, 2, 3]
        mock_tqdm_instance = MagicMock()
        mock_tqdm_instance.__iter__ = MagicMock(return_value=iter(test_data))
        mock_tqdm.return_value = mock_tqdm_instance

        list(iter_progress(test_data, total=10))

        mock_tqdm.assert_called_once_with(test_data, total=10, desc=None)

    @patch('circleseeker.utils.progress.tqdm')
    def test_iter_progress_with_desc_parameter(self, mock_tqdm):
        """Test iter_progress with description parameter."""
        test_data = [1, 2, 3]
        mock_tqdm_instance = MagicMock()
        mock_tqdm_instance.__iter__ = MagicMock(return_value=iter(test_data))
        mock_tqdm.return_value = mock_tqdm_instance

        list(iter_progress(test_data, desc="Processing items"))

        mock_tqdm.assert_called_once_with(test_data, total=None, desc="Processing items")

    @patch('circleseeker.utils.progress.tqdm')
    def test_iter_progress_with_all_parameters(self, mock_tqdm):
        """Test iter_progress with all parameters."""
        test_data = [1, 2, 3]
        mock_tqdm_instance = MagicMock()
        mock_tqdm_instance.__iter__ = MagicMock(return_value=iter(test_data))
        mock_tqdm.return_value = mock_tqdm_instance

        list(iter_progress(test_data, total=5, desc="Custom description", enabled=True))

        mock_tqdm.assert_called_once_with(test_data, total=5, desc="Custom description")

    def test_iter_progress_preserves_order(self):
        """Test that iter_progress preserves the order of items."""
        test_data = [3, 1, 4, 1, 5, 9, 2, 6]
        result = list(iter_progress(test_data))
        assert result == test_data

    def test_iter_progress_lazy_evaluation(self):
        """Test that iter_progress works with lazy evaluation."""
        def generate_numbers():
            """Generator function for testing."""
            for i in range(5):
                yield i * 2

        result = list(iter_progress(generate_numbers()))
        expected = [0, 2, 4, 6, 8]
        assert result == expected

    @patch('circleseeker.utils.progress.tqdm')
    def test_iter_progress_multiple_iterations(self, mock_tqdm):
        """Test iter_progress with multiple separate iterations."""
        test_data1 = [1, 2]
        test_data2 = [3, 4]

        mock_tqdm_instance = MagicMock()
        mock_tqdm_instance.__iter__ = MagicMock(side_effect=[iter(test_data1), iter(test_data2)])
        mock_tqdm.return_value = mock_tqdm_instance

        result1 = list(iter_progress(test_data1))
        result2 = list(iter_progress(test_data2))

        assert result1 == test_data1
        assert result2 == test_data2
        assert mock_tqdm.call_count == 2

    def test_iter_progress_type_annotations(self):
        """Test that iter_progress works with type annotations (implicit test)."""
        # This test verifies that the TypeVar T works correctly
        int_data: List[int] = [1, 2, 3]
        str_data: List[str] = ['a', 'b', 'c']

        int_result = list(iter_progress(int_data))
        str_result = list(iter_progress(str_data))

        assert int_result == int_data
        assert str_result == str_data
        assert all(isinstance(x, int) for x in int_result)
        assert all(isinstance(x, str) for x in str_result)

    @patch('circleseeker.utils.progress.tqdm')
    def test_iter_progress_exception_handling(self, mock_tqdm):
        """Test iter_progress exception handling during tqdm creation."""
        test_data = [1, 2, 3]

        # Test various exceptions that might occur
        for exception in [ImportError(), AttributeError(), ValueError(), Exception()]:
            mock_tqdm.side_effect = exception
            result = list(iter_progress(test_data))
            assert result == test_data

    def test_iter_progress_enabled_parameter_types(self):
        """Test iter_progress enabled parameter with different boolean-like values."""
        test_data = [1, 2, 3]

        # Test with explicit True/False
        assert list(iter_progress(test_data, enabled=True)) == test_data
        assert list(iter_progress(test_data, enabled=False)) == test_data

        # Test with truthy/falsy values
        assert list(iter_progress(test_data, enabled=1)) == test_data
        assert list(iter_progress(test_data, enabled=0)) == test_data
        assert list(iter_progress(test_data, enabled="")) == test_data
        assert list(iter_progress(test_data, enabled="non-empty")) == test_data


class TestIterProgressIntegration:
    """Integration tests for iter_progress."""

    def test_iter_progress_real_world_usage(self):
        """Test iter_progress in a realistic usage scenario."""
        # Simulate processing a list of files
        files = [f"file_{i}.txt" for i in range(10)]
        processed_files = []

        for file in iter_progress(files, desc="Processing files"):
            # Simulate some processing
            processed_files.append(f"processed_{file}")

        assert len(processed_files) == 10
        assert processed_files[0] == "processed_file_0.txt"
        assert processed_files[-1] == "processed_file_9.txt"

    def test_iter_progress_with_computation(self):
        """Test iter_progress with actual computation."""
        numbers = list(range(100))
        squares = []

        for num in iter_progress(numbers, desc="Computing squares"):
            squares.append(num ** 2)

        assert len(squares) == 100
        assert squares[0] == 0
        assert squares[10] == 100
        assert squares[-1] == 99 ** 2

    def test_iter_progress_nested_usage(self):
        """Test iter_progress in nested scenarios."""
        outer_data = [range(5) for _ in range(3)]
        results = []

        for inner_range in iter_progress(outer_data, desc="Outer loop"):
            inner_results = []
            for item in iter_progress(inner_range, desc="Inner loop", enabled=False):
                inner_results.append(item * 2)
            results.append(inner_results)

        assert len(results) == 3
        assert all(len(inner) == 5 for inner in results)
        assert results[0] == [0, 2, 4, 6, 8]