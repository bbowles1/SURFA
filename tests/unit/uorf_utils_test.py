from surfa.uorf_utils import check_identity, fasta_codon_search
import unittest

class TestCheckIdentity(unittest.TestCase):

    def test_check_identity(self):
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


class TestFastaCodonSearch(unittest.TestCase):
    """Test suite for fasta_codon_search function."""

    # Basic functionality tests

    def test_frame_zero_basic(self):
        """Test frame 0 with simple sequence."""
        RNA = "AUGUAA"
        result = fasta_codon_search(RNA, 0)
        self.assertEqual(result, ["AUG", "UAA"])

    def test_frame_one_basic(self):
        """Test frame 1 with simple sequence."""
        RNA = "AUGUAA"
        result = fasta_codon_search(RNA, 1)
        self.assertEqual(result, ["UGU"])  # Note: "AA" is incomplete codon

    def test_frame_two_basic(self):
        """Test frame 2 with simple sequence."""
        RNA = "AUGUAA"
        result = fasta_codon_search(RNA, 2)
        self.assertEqual(result, ["GUA"])  # Note: "GU" is incomplete codon

    # All three reading frames comparison

    def test_all_frames_same_sequence(self):
        """Test all three reading frames produce different results."""
        RNA = "AUGUAGUAA"
        frame0 = fasta_codon_search(RNA, 0)
        frame1 = fasta_codon_search(RNA, 1)
        frame2 = fasta_codon_search(RNA, 2)
        
        self.assertEqual(frame0, ["AUG", "UAG", "UAA"])
        self.assertEqual(frame1, ["UGU", "AGU"], msg="Frame 1 mismatch")
        self.assertEqual(frame2, ["GUA", "GUA"], msg="Frame 2 mismatch")

    # Edge cases

    def test_empty_sequence(self):
        """Test empty RNA sequence."""
        result = fasta_codon_search("", 0)
        self.assertEqual(result, [])

    def test_single_nucleotide(self):
        """Test sequence with only one nucleotide."""
        result = fasta_codon_search("A", 0)
        self.assertEqual(result, [])

    def test_two_nucleotides(self):
        """Test sequence with two nucleotides."""
        result = fasta_codon_search("AU", 0)
        self.assertEqual(result, [])

    def test_three_nucleotides_exact_codon(self):
        """Test sequence with exactly one codon."""
        result = fasta_codon_search("AUG", 0)
        self.assertEqual(result, ["AUG"])

    def test_incomplete_codon_at_end(self):
        """Test sequence with incomplete codon at the end."""
        RNA = "AUGU"  # 4 nucleotides = 1 complete + 1 extra
        result = fasta_codon_search(RNA, 0)
        self.assertEqual(result, ["AUG"])

    def test_five_nucleotides(self):
        """Test 5 nucleotides produces 1 codon with 2 left over."""
        RNA = "AUGUU"
        result = fasta_codon_search(RNA, 0)
        self.assertEqual(result, ["AUG"])

    # Realistic biological sequences

    def test_start_stop_codons(self):
        """Test with common start and stop codons."""
        RNA = "AUGUUUGAUCAGUAA"  # Start, Leu, Stop, Ser, Stop
        result = fasta_codon_search(RNA, 0)
        expected = ["AUG", "UUU", "GAU", "CAG", "UAA"]
        # Note: Only complete triplets are returned
        self.assertEqual(len(result), 5)
        self.assertEqual(result[0], "AUG")
        self.assertEqual(result[1], "UUU")
        self.assertEqual(result[2], "GAU")

    def test_longer_sequence(self):
        """Test with longer realistic sequence."""
        RNA = "AUGUUUCUGGAGUUUUAAGGGCCC"
        result = fasta_codon_search(RNA, 0)
        expected = ["AUG", "UUU", "CUG", "GAG", "UUU", "UAA", "GGG", "CCC"]
        # Should return 5 complete codons
        self.assertEqual(len(result), 8)
        self.assertEqual(result[0], "AUG")  # start codon
        self.assertEqual(result[-1], "CCC")  # end codon

    # Case sensitivity

    def test_uppercase_sequence(self):
        """Test uppercase nucleotides."""
        result = fasta_codon_search("AUGUAA", 0)
        self.assertEqual(result, ["AUG", "UAA"])

    def test_lowercase_sequence(self):
        """Test lowercase nucleotides."""
        result = fasta_codon_search("auguaa", 0)
        self.assertEqual(result, ["aug", "uaa"])

    # Invalid inputs

    def test_negative_frame(self):
        """Test negative frame value (should work like modulo)."""
        RNA = "AUGUAA"
        # should raise value error
        with self.assertRaises(ValueError):
            fasta_codon_search(RNA, -1)

    def test_frame_larger_than_sequence(self):
        """Test frame larger than sequence length."""

        # should raise value error
        with self.assertRaises(ValueError):
            fasta_codon_search("AUG", 10)
        

    def test_frame_three(self):
        """Test frame 3."""
        RNA = "AUGUAA"
        with self.assertRaises(ValueError):
            fasta_codon_search(RNA, 3)

    # Validation tests

    def test_codon_length_validation(self):
        """Verify all returned codons (except possibly last) have length 3."""
        RNA = "AUGUUUCUGGAGUUUUAAGGGCCC"
        result = fasta_codon_search(RNA, 0)
        
        for i, codon in enumerate(result[:-1]):
            self.assertEqual(len(codon), 3, f"Codon {i} should have length 3")
        
        # Last codon can be shorter
        self.assertLessEqual(len(result[-1]), 3)


if __name__ == '__main__':
    unittest.main()
