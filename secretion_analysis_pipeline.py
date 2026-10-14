#!/usr/bin/env python3
"""
Secretion System Analysis Pipeline
Complete workflow for identifying and analyzing secretion systems
Authors: [Your Name]
Date: 2024
"""

import os
import sys
import argparse
import json
import subprocess
from collections import defaultdict
import logging
from Bio import SeqIO
from Bio.Blast import NCBIXML
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ============================================================
# SECRETION SYSTEM DATABASE
# ============================================================

SECRETION_SYSTEMS = {
    'T1SS': {
        'name': 'Type 1 Secretion System',
        'description': 'ABC transporter-mediated secretion',
        'proteins': ['TolC', 'MacB', 'LptB', 'MsbA', 'BtuC', 'BtuD', 'BtuF', 
                    'PstS', 'PstA', 'PstC', 'ZnuA', 'ZnuB', 'ZnuC']
    },
    'T2SS': {
        'name': 'Type 2 Secretion System',
        'description': 'Pseudopili-mediated secretion',
        'proteins': ['GspD', 'GspC', 'GspE', 'GspF', 'GspG', 'GspH', 'GspI',
                    'GspJ', 'GspK', 'GspL', 'GspM', 'GspN', 'GspO', 'GspP',
                    'GspQ', 'GspR', 'GspS', 'HofQ', 'HofB', 'HofC', 'HofD']
    },
    'T3SS': {
        'name': 'Type 3 Secretion System',
        'description': 'Injectisome-mediated effector secretion',
        'proteins': ['SsaV', 'SsaT', 'SsaS', 'SsaR', 'SsaQ', 'SsaP', 'SsaO',
                    'SsaN', 'SsaM', 'SsaL', 'SsaK', 'SsaJ', 'SsaI', 'SsaH',
                    'SsaG', 'SsaF', 'SsaE', 'SsaD', 'SsaC', 'SsaB', 'SsaA',
                    'PrgH', 'PrgI', 'PrgJ', 'PrgK', 'SpaO', 'SpaP', 'SpaQ',
                    'SpaR', 'SpaS', 'InvA', 'InvE', 'InvG']
    },
    'T4SS': {
        'name': 'Type 4 Secretion System',
        'description': 'Conjugation-mediated DNA/protein transfer',
        'proteins': ['TraA', 'TraB', 'TraC', 'TraD', 'TraE', 'TraF', 'TraG',
                    'TraH', 'TraI', 'TraJ', 'TraK', 'TraL', 'TraM', 'TraN',
                    'TraO', 'TraP', 'TraQ', 'TraR', 'TraS', 'TraT', 'TraU',
                    'TraV', 'TraW', 'TraX', 'TraY', 'TraZ', 'VirK', 'TrbI']
    },
    'T5SS': {
        'name': 'Type 5 Secretion System',
        'description': 'Autotransporter-mediated secretion',
        'proteins': ['UspA', 'NhhA', 'SadA', 'Invasin', 'SapB', 'UAT1', 'UAT2',
                    'UAT3', 'UAT4', 'UAT5', 'UAT6', 'UAT7', 'UAT8', 'UAT9',
                    'UAT10', 'UAT11', 'AIDA-I', 'Ag43', 'Hbp', 'IgA1', 'EspP']
    },
    'T6SS': {
        'name': 'Type 6 Secretion System',
        'description': 'Contractile bacteriophage-like secretion',
        'proteins': ['Hcp', 'VgrG', 'ImpA', 'ImpH', 'ImpG', 'ImpF', 'ImpE',
                    'ImpJ', 'ImpK', 'ImpL', 'ImpM', 'VasG', 'VasD', 'VasA',
                    'VasB', 'VasC', 'VasE', 'VasF', 'VasH', 'VasI', 'VasJ',
                    'VasK', 'VasL', 'VasM', 'VasN', 'VasO', 'VasP', 'VasQ',
                    'VasR', 'VasS', 'VasT', 'VasU', 'VasV', 'TssA', 'TssB',
                    'TssC', 'TssD', 'TssE', 'TssF', 'TssG', 'TssH', 'TssI']
    },
    'T8SS': {
        'name': 'Type 8 Secretion System',
        'description': 'Curli/amyloid secretion',
        'proteins': ['CsgA', 'CsgB', 'CsgC', 'CsgD', 'CsgE', 'CsgF', 'CsgG',
                    'Crl', 'MlrA']
    },
    'T9SS': {
        'name': 'Type 9 Secretion System',
        'description': 'Gliding motility-associated secretion',
        'proteins': ['PorV', 'PorU', 'PorP', 'PorQ', 'PorN', 'PorO', 'PorR',
                    'PorS', 'PorT', 'PorW', 'PorX', 'PorY', 'PorZ']
    },
    'FlgSS': {
        'name': 'Flagellar Secretion System',
        'description': 'Flagellar assembly and motility',
        'proteins': ['FlgA', 'FlgB', 'FlgC', 'FlgD', 'FlgE', 'FlgF', 'FlgG',
                    'FlgH', 'FlgI', 'FlgJ', 'FlgK', 'FlgL', 'FlgM', 'FlgN',
                    'FliA', 'FliB', 'FliC', 'FliD', 'FliE', 'FliF', 'FliG',
                    'FliH', 'FliI', 'FliJ', 'FliK', 'FliL', 'FliM', 'FliN',
                    'FliO', 'FliP', 'FliQ', 'FliR', 'FliS', 'FliT', 'FliU',
                    'FliV', 'FliW', 'FliX', 'FliY', 'FliZ']
    },
    'T4PSS': {
        'name': 'Type IV Pili System',
        'description': 'Pilus-mediated adhesion and motility',
        'proteins': ['PilA', 'PilB', 'PilC', 'PilD', 'PilE', 'PilF', 'PilG',
                    'PilH', 'PilI', 'PilJ', 'PilK', 'PilL', 'PilM', 'PilN',
                    'PilO', 'PilP', 'PilQ', 'PilR', 'PilS', 'PilT', 'PilU',
                    'PilV', 'PilW', 'PilX', 'PilY', 'PilZ', 'FimA', 'PpdD',
                    'PpdC', 'PpdA', 'YrfA', 'YrfB', 'YrfC', 'YrfD']
    },
    'SecSS': {
        'name': 'Sec Pathway',
        'description': 'General secretory pathway',
        'proteins': ['SecA', 'SecB', 'SecC', 'SecD', 'SecE', 'SecF', 'SecG',
                    'SecH', 'SecI', 'SecJ', 'SecK', 'SecL', 'SecM', 'SecN',
                    'SecO', 'SecP', 'SecQ', 'SecR', 'SecS', 'SecT', 'SecU',
                    'SecV', 'SecW', 'SecX', 'SecY', 'SecZ', 'MsyB', 'Syd']
    },
    'TatSS': {
        'name': 'Tat Pathway',
        'description': 'Twin-arginine translocation',
        'proteins': ['TatA', 'TatB', 'TatC', 'TatD', 'TatE', 'TatF', 'SufI',
                    'TorA', 'TorD', 'TorE', 'TorF', 'TorG', 'TorH', 'TorI', 'TorJ']
    }
}

# Protein-protein interactions
INTERACTIONS = {
    'T1SS': [
        ('TolC', 'MacB'), ('TolC', 'LptB'), ('TolC', 'MsbA'),
        ('MacB', 'BtuC'), ('BtuC', 'BtuD'), ('BtuC', 'BtuF'),
        ('ZnuA', 'ZnuB'), ('ZnuA', 'ZnuC')
    ],
    'T2SS': [
        ('GspD', 'HofB'), ('GspD', 'HofQ')
    ],
    'T3SS': [
        ('SsaV', 'SsaT'), ('SsaV', 'SsaS'), ('SsaV', 'SsaR'),
        ('PrgH', 'PrgI'), ('PrgH', 'PrgJ'), ('PrgH', 'PrgK'),
        ('InvA', 'InvG'), ('InvA', 'InvE')
    ],
    'T4SS': [
        ('TraA', 'TraB'), ('TraA', 'TraC'), ('TraA', 'TraD'),
        ('TraG', 'TraH'), ('TraG', 'TraI'), ('TraG', 'TraJ')
    ],
    'T5SS': [
        ('UspA', 'NhhA'), ('SadA', 'SapB'), ('Invasin', 'UspA')
    ],
    'T6SS': [
        ('Hcp', 'VgrG'), ('Hcp', 'ImpA'), ('VgrG', 'ImpA'),
        ('ImpA', 'ImpH'), ('VasD', 'VasG')
    ],
    'T8SS': [
        ('CsgA', 'CsgB'), ('CsgB', 'CsgC'), ('CsgC', 'CsgD'),
        ('CsgD', 'CsgE'), ('CsgE', 'CsgF'), ('CsgF', 'CsgG')
    ],
    'FlgSS': [
        ('FlgE', 'FlgF'), ('FlgE', 'FlgG'), ('FliC', 'FliD'),
        ('FliF', 'FliG'), ('FliG', 'FliH'), ('FliM', 'FliN')
    ],
    'T4PSS': [
        ('PilB', 'PilC'), ('PilB', 'PilQ'), ('PilT', 'PilU'),
        ('FimA', 'PilB')
    ],
    'SecSS': [
        ('SecA', 'SecB'), ('SecA', 'SecD'), ('SecD', 'SecE'),
        ('SecE', 'SecF'), ('SecF', 'SecG'), ('SecD', 'SecY')
    ],
    'TatSS': [
        ('TatA', 'TatB'), ('TatA', 'TatC'), ('TatB', 'TatC'),
        ('TatC', 'TatD'), ('TatD', 'TatE')
    ]
}

# ============================================================
# MAIN CLASS
# ============================================================

class SecretionSystemAnalyzer:
    """Main class for secretion system analysis"""
    
    def __init__(self, genome_file, output_dir='results'):
        self.genome_file = genome_file
        self.output_dir = output_dir
        self.genome_name = os.path.basename(genome_file).replace('.fna', '').replace('.fasta', '')
        self.genome = None
        self.contigs = []
        self.tblastn_results = {}
        self.found_proteins = defaultdict(list)
        self.all_hits = {}
        
        # Create output directories
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(f"{output_dir}/reports", exist_ok=True)
        os.makedirs(f"{output_dir}/networks", exist_ok=True)
        os.makedirs(f"{output_dir}/synteny", exist_ok=True)
        os.makedirs(f"{output_dir}/blast", exist_ok=True)
    
    def load_genome(self):
        """Load genome sequences from FASTA file"""
        try:
            self.genome = list(SeqIO.parse(self.genome_file, "fasta"))
            self.contigs = [record.id for record in self.genome]
            logger.info(f"Loaded {len(self.genome)} contigs from {self.genome_file}")
            return self.genome
        except Exception as e:
            logger.error(f"Failed to load genome: {e}")
            raise
    
    def run_tblastn(self, protein_sequences):
        """Run tBLASTn to identify secretion system proteins"""
        import subprocess
        
        logger.info(f"Running tBLASTn for {len(protein_sequences)} proteins...")
        
        # Create query file
        query_file = f"{self.output_dir}/blast/queries.faa"
        with open(query_file, 'w') as f:
            for protein, seq in protein_sequences.items():
                f.write(f">{protein}\n{seq}\n")
        
        # Create BLAST database
        db_name = f"{self.output_dir}/blast/genome_db"
        cmd = f"makeblastdb -in {self.genome_file} -dbtype nucl -out {db_name}"
        subprocess.run(cmd, shell=True, check=False)
        
        # Run tBLASTn
        output_file = f"{self.output_dir}/blast/tblastn_results.txt"
        cmd = f"tblastn -query {query_file} -db {db_name} -evalue 0.01 -outfmt 6 -out {output_file} -max_target_seqs 10 -num_threads 4"
        
        try:
            subprocess.run(cmd, shell=True, check=True)
            logger.info(f"tBLASTn completed successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"tBLASTn failed: {e}")
            return {}
        
        # Parse results
        results = defaultdict(list)
        if os.path.exists(output_file):
            with open(output_file) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        protein = parts[0]
                        hit = {
                            'contig': parts[1],
                            'identity': float(parts[2]),
                            'align_length': int(parts[3]),
                            'evalue': float(parts[4]),
                            'qstart': int(parts[5]),
                            'qend': int(parts[6]),
                            'sstart': int(parts[7]),
                            'send': int(parts[8])
                        }
                        results[protein].append(hit)
        
        logger.info(f"Found hits for {len(results)} proteins")
        return results
    
    def analyze(self, protein_sequences):
        """Run complete analysis"""
        logger.info(f"Starting analysis of {self.genome_name}")
        self.load_genome()
        
        # Run tBLASTn
        self.tblastn_results = self.run_tblastn(protein_sequences)
        
        # Map to secretion systems
        self.found_proteins = defaultdict(list)
        self.all_hits = {}
        
        for system, info in SECRETION_SYSTEMS.items():
            for protein in info['proteins']:
                if protein in self.tblastn_results:
                    self.found_proteins[system].append(protein)
                    self.all_hits[protein] = self.tblastn_results[protein]
                    logger.info(f"  Found {protein} in {system}")
        
        # Generate reports
        self.generate_report()
        self.create_network()
        self.create_synteny()
        
        return dict(self.found_proteins)
    
    def generate_report(self):
        """Generate comprehensive report"""
        report = []
        report.append("="*80)
        report.append(f"SECRETION SYSTEM ANALYSIS REPORT")
        report.append(f"Genome: {self.genome_name}")
        report.append("="*80)
        
        total = sum(len(p) for p in self.found_proteins.values())
        report.append(f"\nTotal proteins identified: {total}")
        report.append(f"Systems detected: {len([s for s in self.found_proteins.keys() if self.found_proteins[s]])}\n")
        
        for system, proteins in sorted(self.found_proteins.items()):
            if proteins:
                report.append(f"\n{system}: {len(proteins)} proteins")
                report.append(f"  {', '.join(proteins[:10])}")
                if len(proteins) > 10:
                    report.append(f"  ... and {len(proteins)-10} more")
        
        # Save report
        with open(f"{self.output_dir}/reports/analysis_report.txt", 'w') as f:
            f.write("\n".join(report))
        
        print("\n".join(report))
    
    def create_network(self):
        """Create protein interaction network"""
        G = nx.Graph()
        
        # Add nodes
        for system, proteins in self.found_proteins.items():
            for protein in proteins:
                G.add_node(protein, system=system)
        
        # Add edges
        for system, interactions in INTERACTIONS.items():
            if system in self.found_proteins:
                for p1, p2 in interactions:
                    if p1 in self.found_proteins[system] and p2 in self.found_proteins[system]:
                        G.add_edge(p1, p2)
        
        # Visualize
        if len(G.nodes()) > 0:
            plt.figure(figsize=(12, 8))
            pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
            
            # Color nodes by system
            colors = []
            for node in G.nodes():
                system = G.nodes[node].get('system', 'Unknown')
                color_map = {
                    'T1SS': '#FF6B6B', 'T2SS': '#4ECDC4', 'T3SS': '#45B7D1',
                    'T4SS': '#96CEB4', 'T5SS': '#FFEAA7', 'T6SS': '#DDA0DD',
                    'T8SS': '#FF9FF3', 'T9SS': '#C9B1FF', 'FlgSS': '#54A0FF',
                    'T4PSS': '#A8E6CF', 'SecSS': '#FFB347', 'TatSS': '#FD79A8'
                }
                colors.append(color_map.get(system, '#DFE6E9'))
            
            nx.draw(G, pos, with_labels=True, node_color=colors,
                   node_size=2000, font_size=8, font_weight='bold',
                   edge_color='gray', width=1.5, alpha=0.7)
            
            plt.title(f"{self.genome_name}\nProtein Interaction Network")
            plt.tight_layout()
            plt.savefig(f"{self.output_dir}/networks/network.png", dpi=300)
            plt.close()
            logger.info(f"Network saved to {self.output_dir}/networks/network.png")
    
    def create_synteny(self):
        """Create synteny plot"""
        # Simplified synteny visualization
        fig, ax = plt.subplots(figsize=(14, 8))
        
        y_pos = 0
        systems_plotted = []
        
        for system, proteins in self.found_proteins.items():
            if proteins:
                y_pos += 1
                systems_plotted.append(system)
                
                # Plot genes
                for i, protein in enumerate(proteins[:15]):
                    x_pos = i * 0.8
                    ax.arrow(x_pos, y_pos, 0.6, 0,
                            head_width=0.15, head_length=0.08,
                            fc='#45B7D1', ec='black', alpha=0.8)
                    ax.text(x_pos + 0.3, y_pos + 0.2, protein,
                           ha='center', va='bottom', fontsize=6, rotation=45)
        
        ax.set_yticks(range(1, len(systems_plotted)+1))
        ax.set_yticklabels(systems_plotted)
        ax.set_xlabel('Gene Order')
        ax.set_title(f"{self.genome_name}\nSynteny of Secretion Systems")
        ax.grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/synteny/synteny.png", dpi=300)
        plt.close()
        logger.info(f"Synteny saved to {self.output_dir}/synteny/synteny.png")

# ============================================================
# MAIN FUNCTION
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Secretion System Analysis Pipeline")
    parser.add_argument('-g', '--genome', required=True, help='Genome FASTA file')
    parser.add_argument('-o', '--output', default='results', help='Output directory')
    parser.add_argument('-p', '--proteins', help='Protein sequences file (optional)')
    
    args = parser.parse_args()
    
    # Load protein sequences (example)
    protein_sequences = {
        'TolC': 'MKKLLPILVLSVAVAVAAAVSLAQPTTIKVKQTPSEQVQTAQNDYYTVYDPAVMTDDAEIEVLTNTF',
        'MacB': 'MSLKRLLLSIAAVLVAAVSLATPVVQSAKDVSAVVKETTTEVVKDVSVK',
        'Hcp': 'MIIKDLTIAKNLPTELKAAINAAITSLKDEISQLTAKFEKAIIAQ',
        'VgrG': 'MKKILFSAVMALSLATALATAAVDVNKVTAVAKTLKETTAKDVAAK',
        'SecA': 'MKKILVFTTLLAIGAAVVTAAVDAVKETTVKEVTAKVEKATVKET',
    }
    
    analyzer = SecretionSystemAnalyzer(args.genome, args.output)
    results = analyzer.analyze(protein_sequences)
    
    print(f"\n✅ Analysis complete! Results in {args.output}")

if __name__ == "__main__":
    main()
