#TERMINAL COMMAND TO RUN THIS CODE:
#python GeneAnalysis.py Zea_mays.B73_RefGen_v4.48.chr.gff3.gz Zea_mays.B73_RefGen_v4.dna.chromosome.1.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.2.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.3.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.4.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.5.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.6.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.7.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.8.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.9.fa.gz Zea_mays.B73_RefGen_v4.dna.chromosome.10.fa.gz promoters.txt zea_mays_genes.txt
#Before running the code please enzure that all the files mentioned in the command are in the same directory

import argparse  # For parsing command-line arguments
from Bio import SeqIO  # Biopython library for parsing sequence files
import gzip  # For reading compressed files
import random  # For generating random samples
import re  # Regular expressions for pattern matching
from typing import Dict, List, Set, Tuple  # Type hints for better code readability and type checking

def parse_arguments():
    """
    Parse command-line arguments for the script.
    
    This function sets up an argument parser to handle input files and optional parameters:
    - gff_file: Compressed GFF3 file with gene annotations
    - genome_files: Compressed FASTA files for chromosomes
    - motifs_file: File containing promoter motifs to search
    - target_genes_file: File with list of target genes to analyze
    - promoter_length: Optional parameter to define promoter region size (default 500 bp)
    - random_samples: Optional parameter to set number of random gene samples (default 5)
    
    Returns:
        Parsed arguments object with all input parameters
    """
    parser = argparse.ArgumentParser(description='Analyze promoter regions in maize genes')
    parser.add_argument('gff_file', help='Compressed GFF3 file containing gene annotations')
    parser.add_argument('genome_files', nargs='+', help='Compressed FASTA files containing genome sequences (chromosomes 1-10)')
    parser.add_argument('motifs_file', help='File containing promoter motifs')
    parser.add_argument('target_genes_file', help='File containing list of target genes')
    parser.add_argument('--promoter_length', type=int, default=500,
                       help='Length of promoter region to analyze (default: 500)')
    parser.add_argument('--random_samples', type=int, default=5,
                       help='Number of random gene samples to analyze (default: 5)')
    return parser.parse_args()

def read_genome(genome_files: List[str]) -> Dict[str, str]:
    """
    Read genome sequences from compressed FASTA files.
    
    This function reads multiple genome files (chromosomes) and creates a dictionary
    mapping chromosome names to their complete sequences. Key features:
    - Uses gzip to read compressed files
    - Uses Biopython's SeqIO for parsing FASTA formats
    - Handles exceptions for file reading
    - Converts sequences to strings for easier manipulation
    
    Args:
        genome_files (List[str]): Paths to compressed genome FASTA files
    
    Returns:
        Dict[str, str]: Dictionary of chromosome names and their sequences
    
    Raises:
        ValueError if no sequences are found in any files
    """
    genome = {}
    for genome_file in genome_files:
        try:
            with gzip.open(genome_file, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    genome[record.id] = str(record.seq)
        except Exception as e:
            print(f"Error reading genome file {genome_file}: {e}")
            raise
    
    if not genome:
        raise ValueError("No sequences found in genome files")
    return genome

def parse_gff(gff_file: str, genome_files: List[str]) -> Dict[str, List[Tuple[str, int, str]]]:
    """
    Parse GFF3 file to extract gene information for specific chromosomes.
    
    Complex function with multiple important steps:
    1. Extract chromosome numbers from genome file names
    2. Skip non-gene entries and comments
    3. Validate and parse GFF3 lines
    4. Extract key gene information:
       - Chromosome
       - Start position
       - Strand (+ or -)
       - Gene ID
    
    Robust error handling includes:
    - Skipping malformed lines
    - Handling missing attributes
    - Only processing specified chromosomes
    
    Args:
        gff_file (str): Path to compressed GFF3 file
        genome_files (List[str]): Genome files to determine valid chromosomes
    
    Returns:
        Dict[str, List[Tuple[str, int, str]]]: 
        Gene IDs mapped to lists of (chromosome, start_position, strand)
    """
    # Extract chromosome numbers from genome file names
    valid_chromosomes = [
        file.split('chromosome.')[1].split('.fa.gz')[0] 
        for file in genome_files
    ]
    
    genes = {}
    try:
        with gzip.open(gff_file, 'rt') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    print(f"Warning: Skipping malformed line: {line.strip()}")
                    continue
                    
                if fields[2] != 'gene':
                    continue
                
                # Only process chromosomes matching the genome files
                if fields[0] not in valid_chromosomes:
                    continue
                
                try:
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    
                    # Parse attributes
                    attributes = {}
                    for attr in fields[8].split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attributes[key] = value
                    
                    gene_id = attributes.get('ID', '').split(':')[-1]
                    
                    if gene_id:
                        if gene_id not in genes:
                            genes[gene_id] = []
                        genes[gene_id].append((chrom, start, strand))
                except Exception as e:
                    print(f"Warning: Error parsing line: {line.strip()}")
                    print(f"Error details: {e}")
                    continue
                    
        if not genes:
            raise ValueError(f"No genes found in {gff_file}")
            
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        raise
        
    print(f"Genes on specified chromosomes: {len(genes)}")
    return genes

def read_target_genes(target_genes_file: str) -> Set[str]:
    """
    Read list of target genes from a file.
    
    Simple function to:
    - Read gene names from a file
    - Convert to a set for fast lookup
    - Handle empty files
    
    Args:
        target_genes_file (str): Path to file with target gene names
    
    Returns:
        Set[str]: Unique gene names
    
    Raises:
        ValueError if no genes are found
    """
    try:
        with open(target_genes_file) as f:
            genes = {line.strip() for line in f if line.strip()}
        if not genes:
            raise ValueError(f"No genes found in {target_genes_file}")
        return genes
    except Exception as e:
        print(f"Error reading target genes file: {e}")
        raise

def read_motifs(motifs_file: str) -> Dict[str, str]:
    """
    Read promoter motifs from a file.
    
    Handles:
    - Reading motif patterns
    - Converting to uppercase
    - Skipping empty lines
    - Providing diagnostic information about read motifs
    
    Args:
        motifs_file (str): Path to file with motif patterns
    
    Returns:
        Dict[str, str]: Motif names mapped to their patterns
    
    Raises:
        ValueError if no motifs are found
    """
    try:
        motifs = {}
        with open(motifs_file) as f:
            for line in f:
                motif = line.strip().upper()  # Convert to uppercase
                if motif:  # Skip empty lines
                   motifs[motif] = motif
        if not motifs:
            raise ValueError(f"No motifs found in {motifs_file}")
        print(f"Successfully read {len(motifs)} motifs:")
        for motif_id, pattern in list(motifs.items())[:5]:  # Show first 5 motifs
            print(f"  {motif_id}: {pattern}")
        if len(motifs) > 5:
            print(f"  ... and {len(motifs)-5} more motifs")
        return motifs
    except Exception as e:
        print(f"Error reading motifs file: {e}")
        raise

def count_motifs(sequence: str, motifs: Dict[str, str]) -> Dict[str, int]:
    """
    Count occurrences of motifs in a given sequence.
    
    Key points:
    - Case-insensitive matching
    - Uses regex to find all motif occurrences
    - Handles potential regex compilation errors
    
    Args:
        sequence (str): DNA sequence to search
        motifs (Dict[str, str]): Dictionary of motif patterns
    
    Returns:
        Dict[str, int]: Motif names mapped to their count in the sequence
    """
    sequence = sequence.upper()  # Convert sequence to uppercase for matching
    counts = {}
    for motif_name, motif_pattern in motifs.items():
        try:
            counts[motif_name] = len(re.findall(motif_pattern, sequence))
        except re.error:
            print(f"Warning: Invalid regex pattern for motif {motif_name}: {motif_pattern}")
            counts[motif_name] = 0
    return counts

def extract_promoter_sequence(genome: Dict[str, str], chrom: str, pos: int, 
                            strand: str, length: int) -> str:
    """
    Extract promoter sequence from genome, handling different strands.
    
    Complex function managing:
    - Extracting sequences on positive and negative strands
    - Handling sequence boundaries
    - Reverse complementing for negative strand
    - Removing long N-string regions (assembly gaps)
    
    Args:
        genome (Dict[str, str]): Full genome sequences
        chrom (str): Chromosome name
        pos (int): Transcription start site position
        strand (str): Gene strand (+ or -)
        length (int): Promoter region length
    
    Returns:
        str: Extracted promoter sequence
    """
    if strand == '+':
        start = max(0, pos - length)
        seq = genome[chrom][start:pos]
    else:  # strand == '-'
        start = pos
        end = min(len(genome[chrom]), pos + length)
        seq = genome[chrom][start:end]
        # Convert to reverse complement
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        seq = ''.join(complement.get(base, base) for base in reversed(seq))
    
    # Check for N-strings (assembly gaps)
    n_position = seq.find('N' * 100)
    if n_position != -1:
        seq = seq[:n_position]
    
    return seq

def analyze_gene_set(genes: Dict[str, List[Tuple[str, int, str]]], 
                    gene_ids: Set[str], genome: Dict[str, str],
                    motifs: Dict[str, str], promoter_length: int) -> Dict[str, int]:
    """
    Comprehensive gene set promoter analysis.
    
    Performs:
    - Gene promoter sequence extraction
    - Motif counting across promoter regions
    - Progress tracking
    
    Args:
        genes (Dict): Mapping of gene IDs to genomic locations
        gene_ids (Set): Set of genes to analyze
        genome (Dict): Full genome sequences
        motifs (Dict): Motif patterns to search
        promoter_length (int): Length of promoter region
    
    Returns:
        Dict[str, int]: Total motif counts across analyzed genes
    """
    total_counts = {motif: 0 for motif in motifs}
    processed_genes = 0
    total_genes = len(gene_ids)
    
    for gene_id in gene_ids:
        if gene_id in genes:
            # Use most upstream TSS if multiple exist
            gene_info = min(genes[gene_id], key=lambda x: x[1])
            chrom, pos, strand = gene_info
            
            promoter_seq = extract_promoter_sequence(genome, chrom, pos, 
                                                   strand, promoter_length)
            
            counts = count_motifs(promoter_seq, motifs)
            for motif, count in counts.items():
                total_counts[motif] += count
            
            processed_genes += 1
            if processed_genes % 100 == 0:
                print(f"Processed {processed_genes}/{total_genes} genes...")
    
    return total_counts

def main():
    """
    Main script execution function.
    
    Orchestrates entire analysis workflow:
    1. Parse command-line arguments
    2. Read genome files
    3. Parse gene annotations
    4. Read target and random gene sets
    5. Analyze motif representation
    6. Write comprehensive results
    
    Error handling and logging throughout the process
    """
    try:
        args = parse_arguments()
        
        print("\nReading input files...")
        print(f"Processing GFF file: {args.gff_file}")
        genes = parse_gff(args.gff_file, args.genome_files)
        print(f"Found {len(genes)} genes")
        
        print(f"\nProcessing genome files: {args.genome_files}")
        genome = read_genome(args.genome_files)
        print(f"Found {len(genome)} chromosomes")
        
        print(f"\nReading target genes: {args.target_genes_file}")
        target_genes = read_target_genes(args.target_genes_file)
        print(f"Found {len(target_genes)} target genes")
        
        print(f"\nReading motifs: {args.motifs_file}")
        motifs = read_motifs(args.motifs_file)
        
        print("\nAnalyzing target genes...")
        target_counts = analyze_gene_set(genes, target_genes, genome, motifs, 
                                       args.promoter_length)
        
        print("\nAnalyzing random gene sets...")
        all_genes = set(genes.keys())
        random_counts = []
        for i in range(args.random_samples):
            print(f"Processing random set {i+1}/{args.random_samples}")
            random_genes = random.sample(list(all_genes), len(target_genes))

            counts = analyze_gene_set(genes, random_genes, genome, motifs,
                                    args.promoter_length)
            random_counts.append(counts)
        
        print("\nWriting results...")
        # Write detailed results
        with open('motif_counts.txt', 'w') as f:
            f.write('Motif Count Analysis:\n\n')
            header = ['Motif', 'Counts_in_Selected_Genes'] + [f'Counts_in_Random_Set_{i+1}' 
                                                 for i in range(args.random_samples)]
            f.write('\t'.join(header) + '\n')
            
            for motif in motifs:
                row = [motif, str(target_counts[motif])]
                row.extend(str(counts[motif]) for counts in random_counts)
                f.write('\t'.join(row) + '\n')
        
        print("\nAnalyzing motif representation...")
        # Statistically assess motif representation
        with open('motif_analysis.txt', 'w') as f:
            f.write('Motif Analysis:\n\n')
            for motif in motifs:
                target_count = target_counts[motif]
                random_mean = sum(counts[motif] for counts in random_counts) / args.random_samples
                
                if target_count > 2 * random_mean:
                    f.write(f'{motif}: over-represented '
                           f'(Target: {target_count}, Random avg: {random_mean:.2f})\n')
                elif target_count < 0.5 * random_mean:
                    f.write(f'{motif}: under-represented '
                           f'(Target: {target_count}, Random avg: {random_mean:.2f})\n')
        
        print("\nAnalysis complete! Results written to motif_counts.txt and motif_analysis.txt")
                
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        raise

if __name__ == '__main__':
    main()