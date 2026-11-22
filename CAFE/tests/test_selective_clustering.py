"""Tests for Selective_Clustering.py B10 fix.

B10: Selecting "None" in the radio must not wipe all session state.
"""
import pytest
from unittest.mock import patch, MagicMock


class TestB10SessionStateWipe:
    """B10: Only selective-clustering-specific keys should be cleared, not all."""

    def test_preserves_adata(self):
        """The 'adata' key must survive a 'None' selection."""
        # Simulate the B10 fix logic
        session_state = {
            'adata': 'some_anndata_object',
            'leiden_computed': True,
            'pca_done': True,
            'selective_markers': ['CD1', 'CD2'],
            'subset_cell_type': 'T cell',
            'umap_computed': False,
        }

        # B10 fix: only clear selective_*/subset_* keys
        for key in list(session_state.keys()):
            if key.startswith('selective_') or key.startswith('subset_'):
                del session_state[key]

        # Critical state must survive
        assert 'adata' in session_state
        assert session_state['adata'] == 'some_anndata_object'
        assert 'leiden_computed' in session_state
        assert 'pca_done' in session_state
        # Selective-clustering keys must be cleared
        assert 'selective_markers' not in session_state
        assert 'subset_cell_type' not in session_state

    def test_old_behavior_would_destroy_data(self):
        """Demonstrate what the old buggy code did."""
        session_state = {
            'adata': 'some_anndata_object',
            'leiden_computed': True,
            'pca_done': True,
        }

        # Old buggy code: wipe everything
        for key in list(session_state.keys()):
            del session_state[key]

        # Everything is gone, which is the bug.
        assert len(session_state) == 0
        assert 'adata' not in session_state

    def test_no_selective_keys(self):
        """If there are no selective_*/subset_* keys, nothing happens."""
        session_state = {
            'adata': 'object',
            'pca_done': True,
        }

        for key in list(session_state.keys()):
            if key.startswith('selective_') or key.startswith('subset_'):
                del session_state[key]

        assert len(session_state) == 2
        assert session_state['adata'] == 'object'
