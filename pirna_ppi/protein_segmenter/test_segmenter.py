#!/usr/bin/env python3
"""Test script for protein segmenter module."""

import sys
from pathlib import Path


def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")
    try:
        from pirna_ppi.protein_segmenter import Domain, Protein, Segment, SegmenterConfig
        from pirna_ppi.protein_segmenter.core import ProteinSegmenter, calculate_inter_element_contacts
        from pirna_ppi.protein_segmenter.parsers import (
            parse_dpam_domains,
            parse_fasta,
            parse_iprscan,
            parse_pae,
            parse_pdb_coordinates,
        )

        print("  All imports successful!")
        return True
    except ImportError as e:
        print(f"  Import error: {e}")
        return False


def test_domain_model():
    """Test Domain model."""
    print("Testing Domain model...")
    from pirna_ppi.protein_segmenter.models import Domain

    # Test single range
    d1 = Domain(name="nD1", ranges=[(1, 100)])
    assert d1.total_length == 100, f"Expected 100, got {d1.total_length}"
    assert d1.start == 1, f"Expected 1, got {d1.start}"
    assert d1.end == 100, f"Expected 100, got {d1.end}"

    # Test multiple ranges
    d2 = Domain(name="nD2", ranges=[(1, 50), (100, 150)])
    # (50-1+1) + (150-100+1) = 50 + 51 = 101
    assert d2.total_length == 101, f"Expected 101, got {d2.total_length}"
    assert d2.start == 1, f"Expected 1, got {d2.start}"
    assert d2.end == 150, f"Expected 150, got {d2.end}"

    print("  Domain model tests passed!")
    return True


def test_range_parser():
    """Test range string parsing."""
    print("Testing range parser...")
    from pirna_ppi.protein_segmenter.parsers.domain_parser import parse_range_string

    # Test simple range
    ranges = parse_range_string("1-100")
    assert ranges == [(1, 100)]

    # Test multiple ranges
    ranges = parse_range_string("1-50,100-150,200-250")
    assert ranges == [(1, 50), (100, 150), (200, 250)]

    print("  Range parser tests passed!")
    return True


def test_merge_criteria():
    """Test element pair merge criteria."""
    print("Testing merge criteria...")
    from pirna_ppi.protein_segmenter.models import Domain, ElementPairMetrics

    d1 = Domain(name="nD1", ranges=[(1, 100)])
    d2 = Domain(name="nD2", ranges=[(150, 250)])

    # Case 1: High contacts, low PAE -> should merge
    metrics1 = ElementPairMetrics(
        element1=d1,
        element2=d2,
        n_contact=50,
        l_min=100,
        avg_pae=5.0,
    )
    # 50/100 + 50/100 = 1.0 > 5.0/8 = 0.625 -> should merge
    assert metrics1.should_merge is True

    # Case 2: Low contacts, high PAE -> should not merge
    metrics2 = ElementPairMetrics(
        element1=d1,
        element2=d2,
        n_contact=5,
        l_min=100,
        avg_pae=20.0,
    )
    # 5/100 + 5/100 = 0.1 < 20.0/8 = 2.5 -> should not merge
    assert metrics2.should_merge is False

    print("  Merge criteria tests passed!")
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Protein Segmenter Module Tests")
    print("=" * 60)

    tests = [
        test_imports,
        test_domain_model,
        test_range_parser,
        test_merge_criteria,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  Test failed with exception: {e}")
            failed += 1

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
