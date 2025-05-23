import re
from typing import List, Tuple, Dict

class RestrictionEnzyme:
    """Class to represent a restriction enzyme with its recognition sequence and cut sites."""
    
    def __init__(self, name: str, recognition_seq: str, cut_pos_5: int, cut_pos_3: int = None):
        """
        Initialize restriction enzyme.
        
        Args:
            name: Name of the enzyme (e.g., 'EcoRI')
            recognition_seq: Recognition sequence (e.g., 'GAATTC')
            cut_pos_5: Cut position on 5' strand (0-based from start of recognition seq)
            cut_pos_3: Cut position on 3' strand (if None, assumes blunt end)
        """
        self.name = name
        self.recognition_seq = recognition_seq.upper()
        self.cut_pos_5 = cut_pos_5
        self.cut_pos_3 = cut_pos_3 if cut_pos_3 is not None else cut_pos_5
        
    def get_overhang_type(self):
        """Determine if enzyme creates sticky (cohesive) or blunt ends."""
        if self.cut_pos_5 == self.cut_pos_3:
            return "blunt"
        elif self.cut_pos_5 > self.cut_pos_3:
            return "5' overhang"
        else:
            return "3' overhang"

class DNADigest:
    """Class to perform restriction enzyme digestion analysis."""
    
    def __init__(self):
        self.custom_enzymes = {}
    
    def add_custom_enzyme(self, name: str, recognition_seq: str, cut_pos_5: int, cut_pos_3: int = None):
        """
        Add a custom restriction enzyme.
        
        Args:
            name: Name of the enzyme
            recognition_seq: Recognition sequence
            cut_pos_5: Cut position on 5' strand (0-based from start of recognition seq)
            cut_pos_3: Cut position on 3' strand (if None, assumes same as cut_pos_5)
        """
        enzyme = RestrictionEnzyme(name, recognition_seq, cut_pos_5, cut_pos_3)
        self.custom_enzymes[name] = enzyme
        return enzyme
    
    def find_restriction_sites(self, dna_seq: str, recognition_seq: str) -> List[int]:
        """
        Find all restriction sites for a given recognition sequence in DNA sequence.
        
        Args:
            dna_seq: DNA sequence string
            recognition_seq: Recognition sequence to search for
            
        Returns:
            List of positions where restriction sites are found
        """
        dna_seq = dna_seq.upper().replace(' ', '').replace('\n', '').replace('\t', '')
        recognition_seq = recognition_seq.upper()
        
        # Validate DNA sequence
        valid_bases = set('ATCG')
        if not all(base in valid_bases for base in dna_seq):
            raise ValueError("DNA sequence contains invalid characters. Use only A, T, C, G.")
        
        # Validate recognition sequence
        if not all(base in valid_bases for base in recognition_seq):
            raise ValueError("Recognition sequence contains invalid characters. Use only A, T, C, G.")
        
        # Find all occurrences of recognition sequence
        sites = []
        for i in range(len(dna_seq) - len(recognition_seq) + 1):
            if dna_seq[i:i+len(recognition_seq)] == recognition_seq:
                sites.append(i)
        
        return sites
    
    def digest_dna(self, dna_seq: str, recognition_seq: str, cut_pos_5: int, cut_pos_3: int = None, enzyme_name: str = "Custom") -> Dict:
        """
        Perform restriction digest and return fragments.
        
        Args:
            dna_seq: DNA sequence string
            recognition_seq: Recognition sequence
            cut_pos_5: Cut position on 5' strand (0-based from start of recognition seq)
            cut_pos_3: Cut position on 3' strand (if None, assumes same as cut_pos_5)
            enzyme_name: Name for the enzyme (optional)
            
        Returns:
            Dictionary containing digest results
        """
        dna_seq = dna_seq.upper().replace(' ', '').replace('\n', '').replace('\t', '')
        sites = self.find_restriction_sites(dna_seq, recognition_seq)
        
        if cut_pos_3 is None:
            cut_pos_3 = cut_pos_5
        
        # Create temporary enzyme object for overhang type calculation
        temp_enzyme = RestrictionEnzyme(enzyme_name, recognition_seq, cut_pos_5, cut_pos_3)
        
        if not sites:
            return {
                'enzyme': enzyme_name,
                'recognition_sequence': recognition_seq,
                'cut_sites': [],
                'fragments': [dna_seq],
                'fragment_lengths': [len(dna_seq)],
                'overhang_type': temp_enzyme.get_overhang_type(),
                'total_sites': 0
            }
        
        # Calculate actual cut positions (using 5' strand cut position)
        cut_positions = []
        for site in sites:
            cut_pos = site + cut_pos_5
            if 0 <= cut_pos <= len(dna_seq):  # Ensure cut position is valid
                cut_positions.append(cut_pos)
        
        # Generate fragments
        fragments = []
        start = 0
        
        for cut_pos in sorted(cut_positions):
            if start < cut_pos:
                fragments.append(dna_seq[start:cut_pos])
            start = cut_pos
        
        # Add final fragment
        if start < len(dna_seq):
            fragments.append(dna_seq[start:])
        
        # Remove empty fragments
        fragments = [f for f in fragments if f]
        
        return {
            'enzyme': enzyme_name,
            'recognition_sequence': recognition_seq,
            'cut_sites': cut_positions,
            'fragments': fragments,
            'fragment_lengths': [len(f) for f in fragments],
            'overhang_type': temp_enzyme.get_overhang_type(),
            'total_sites': len(sites),
            'recognition_sites': sites
        }
    
    def multiple_digest(self, dna_seq: str, enzymes_data: List[Tuple]) -> Dict:
        """
        Perform digest with multiple enzymes.
        
        Args:
            dna_seq: DNA sequence string
            enzymes_data: List of tuples (name, recognition_seq, cut_pos_5, cut_pos_3)
            
        Returns:
            Dictionary containing multiple digest results
        """
        dna_seq = dna_seq.upper().replace(' ', '').replace('\n', '').replace('\t', '')
        
        all_cuts = []
        enzyme_info = []
        
        for enzyme_data in enzymes_data:
            name, recognition_seq, cut_pos_5 = enzyme_data[:3]
            cut_pos_3 = enzyme_data[3] if len(enzyme_data) > 3 else cut_pos_5
            
            sites = self.find_restriction_sites(dna_seq, recognition_seq)
            cut_positions = [site + cut_pos_5 for site in sites if 0 <= site + cut_pos_5 <= len(dna_seq)]
            
            all_cuts.extend(cut_positions)
            enzyme_info.append({
                'name': name,
                'recognition_seq': recognition_seq,
                'sites': len(sites),
                'cut_positions': cut_positions
            })
        
        # Sort all cut positions
        all_cuts = sorted(set(all_cuts))
        
        # Generate fragments
        fragments = []
        start = 0
        
        for cut_pos in all_cuts:
            if start < cut_pos:
                fragments.append(dna_seq[start:cut_pos])
            start = cut_pos
        
        # Add final fragment
        if start < len(dna_seq):
            fragments.append(dna_seq[start:])
        
        return {
            'enzymes': enzyme_info,
            'total_cuts': len(all_cuts),
            'cut_positions': all_cuts,
            'fragments': fragments,
            'fragment_lengths': [len(f) for f in fragments]
        }

def validate_cut_position(recognition_seq: str, cut_pos: int) -> bool:
    """Validate that cut position is within the recognition sequence length."""
    return 0 <= cut_pos <= len(recognition_seq)

def get_user_input():
    """Get DNA sequence and enzyme information from user."""
    print("=== Restriction Enzyme Digestion Analysis ===\n")
    
    # Get DNA sequence
    print("Enter your DNA sequence:")
    print("(You can paste multiple lines, press Enter twice when done)")
    dna_lines = []
    while True:
        line = input().strip()
        if line == "" and dna_lines:
            break
        if line:
            dna_lines.append(line)
    
    dna_sequence = ''.join(dna_lines).upper().replace(' ', '').replace('\n', '').replace('\t', '')
    
    if not dna_sequence:
        print("No DNA sequence provided!")
        return None, None
    
    print(f"\nDNA sequence received ({len(dna_sequence)} bp)")
    if len(dna_sequence) <= 100:
        print(f"Sequence: {dna_sequence}")
    else:
        print(f"Sequence: {dna_sequence[:50]}...{dna_sequence[-50:]}")
    
    # Get enzyme information
    enzymes = []
    
    while True:
        print(f"\n--- Enzyme {len(enzymes) + 1} ---")
        
        # Get enzyme name
        name = input("Enter enzyme name (or press Enter to finish): ").strip()
        if not name:
            break
        
        # Get recognition sequence
        recognition_seq = input("Enter recognition sequence: ").strip().upper()
        if not recognition_seq:
            print("Recognition sequence cannot be empty!")
            continue
        
        # Validate recognition sequence
        if not all(base in 'ATCG' for base in recognition_seq):
            print("Invalid recognition sequence! Use only A, T, C, G.")
            continue
        
        # Get cut position
        while True:
            try:
                cut_input = input(f"Enter cut position on 5' strand (0-{len(recognition_seq)}): ")
                cut_pos_5 = int(cut_input)
                if not validate_cut_position(recognition_seq, cut_pos_5):
                    print(f"Cut position must be between 0 and {len(recognition_seq)}")
                    continue
                break
            except ValueError:
                print("Please enter a valid number!")
        
        # Ask about 3' cut position
        cut_pos_3 = None
        three_prime = input("Different cut position on 3' strand? (y/n): ").strip().lower()
        if three_prime == 'y':
            while True:
                try:
                    cut_input = input(f"Enter cut position on 3' strand (0-{len(recognition_seq)}): ")
                    cut_pos_3 = int(cut_input)
                    if not validate_cut_position(recognition_seq, cut_pos_3):
                        print(f"Cut position must be between 0 and {len(recognition_seq)}")
                        continue
                    break
                except ValueError:
                    print("Please enter a valid number!")
        
        enzymes.append((name, recognition_seq, cut_pos_5, cut_pos_3))
        
        # Ask if user wants to add more enzymes
        more = input("Add another enzyme? (y/n): ").strip().lower()
        if more != 'y':
            break
    
    return dna_sequence, enzymes

def print_digest_results(results: Dict):
    """Pretty print digest results."""
    print(f"\n=== Digest Results: {results['enzyme']} ===")
    print(f"Recognition Sequence: {results['recognition_sequence']}")
    print(f"Overhang Type: {results['overhang_type']}")
    print(f"Recognition Sites Found: {results['total_sites']}")
    
    if results['recognition_sites']:
        print(f"Recognition Site Positions: {results['recognition_sites']}")
        print(f"Cut Positions: {results['cut_sites']}")
    
    print(f"\nFragments Generated: {len(results['fragments'])}")
    for i, (fragment, length) in enumerate(zip(results['fragments'], results['fragment_lengths']), 1):
        print(f"  Fragment {i}: {length} bp")
        if length <= 60:  # Show sequence for shorter fragments
            print(f"    Sequence: {fragment}")
        else:
            print(f"    Sequence: {fragment[:30]}...{fragment[-30:]}")

def print_multiple_digest_results(results: Dict):
    """Pretty print multiple digest results."""
    print(f"\n=== Multiple Enzyme Digest Results ===")
    
    for enzyme in results['enzymes']:
        print(f"{enzyme['name']}: {enzyme['sites']} sites ({enzyme['recognition_seq']})")
    
    print(f"\nTotal Cut Positions: {len(results['cut_positions'])}")
    if results['cut_positions']:
        print(f"Cut Positions: {results['cut_positions']}")
    
    print(f"\nFragments Generated: {len(results['fragments'])}")
    for i, (fragment, length) in enumerate(zip(results['fragments'], results['fragment_lengths']), 1):
        print(f"  Fragment {i}: {length} bp")
        if length <= 60:
            print(f"    Sequence: {fragment}")
        else:
            print(f"    Sequence: {fragment[:30]}...{fragment[-30:]}")

def main():
    """Main function to run the interactive restriction enzyme analysis."""
    
    digest_analyzer = DNADigest()
    
    try:
        dna_sequence, enzymes_data = get_user_input()
        
        if not dna_sequence or not enzymes_data:
            print("No valid input provided. Exiting.")
            return
        
        if len(enzymes_data) == 1:
            # Single enzyme digest
            name, recognition_seq, cut_pos_5, cut_pos_3 = enzymes_data[0]
            results = digest_analyzer.digest_dna(dna_sequence, recognition_seq, cut_pos_5, cut_pos_3, name)
            print_digest_results(results)
        
        else:
            # Multiple enzyme digest
            results = digest_analyzer.multiple_digest(dna_sequence, enzymes_data)
            print_multiple_digest_results(results)
            
            # Also show individual enzyme results
            print(f"\n=== Individual Enzyme Results ===")
            for enzyme_data in enzymes_data:
                name, recognition_seq, cut_pos_5 = enzyme_data[:3]
                cut_pos_3 = enzyme_data[3] if len(enzyme_data) > 3 else cut_pos_5
                individual_result = digest_analyzer.digest_dna(dna_sequence, recognition_seq, cut_pos_5, cut_pos_3, name)
                print_digest_results(individual_result)
    
    except Exception as e:
        print(f"Error: {e}")
        print("Please check your input and try again.")

if __name__ == "__main__":
    main()