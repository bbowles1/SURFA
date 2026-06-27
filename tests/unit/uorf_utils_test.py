import pytest
import sys
import os

def test_check_identity():
    """
    Parametrized test for check_identity function.
    
    Covers logic for:
    1. Positive strand: 5' is upstream (lower coordinate), 3' is downstream (higher coordinate).
    2. Negative strand: 5' is downstream (higher coordinate), 3' is upstream (lower coordinate).
    """
    
    test_data = [
        # (region_start, strand, CDS_start, expected_result, description)
        
        # Positive Strand Cases
        (100, "+", 200, 5, "Positive strand, UTR upstream of CDS (5' UTR)"),
        (300, "+", 200, 3, "Positive strand, UTR downstream of CDS (3' UTR)"),
        
        # Negative Strand Cases
        (300, "-", 200, 5, "Negative strand, UTR downstream of CDS (5' UTR)"),
        (100, "-", 200, 3, "Negative strand, UTR upstream of CDS (3' UTR)"),
        
        # Edge cases for strict inequality
        (200, "+", 200, None, "Positive strand, Region equals CDS start"),
        (200, "-", 200, None, "Negative strand, Region equals CDS start"),
    ]

    for region_start, strand, cds_start, expected, description in test_data:
        result = check_identity(region_start, strand, cds_start)
        assert result == expected, f"Failed for case: {description}. Expected {expected}, got {result}"

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
