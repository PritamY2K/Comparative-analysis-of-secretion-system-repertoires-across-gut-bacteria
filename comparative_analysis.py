#!/usr/bin/env python3
"""
Comparative analysis of secretion systems across multiple genomes
"""

import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict

class ComparativeAnalyzer:
    def __init__(self, results_dir='results'):
        self.results_dir = results_dir
        self.genomes = []
        self.data = {}
        self.load_results()
    
    def load_results(self):
        """Load results from all genome analyses"""
        for genome_dir in os.listdir(self.results_dir):
            report_file = f"{self.results_dir}/{genome_dir}/reports/results_summary.json"
            if os.path.exists(report_file):
                with open(report_file) as f:
                    self.data[genome_dir] = json.load(f)
                    self.genomes.append(genome_dir)
        
        print(f"Loaded results for {len(self.genomes)} genomes")
    
    def create_comparison_table(self):
        """Create comparison table"""
        systems = ['T1SS', 'T2SS', 'T3SS', 'T4SS', 'T5SS', 'T6SS', 
                  'T8SS', 'T9SS', 'FlgSS', 'T4PSS', 'SecSS', 'TatSS']
        
        data = {'System': systems}
        for genome in self.genomes:
            counts = []
            for system in systems:
                proteins = self.data[genome].get('systems_present', {}).get(system, [])
                counts.append(len(proteins))
            data[genome] = counts
        
        df = pd.DataFrame(data)
        df.set_index('System', inplace=True)
        return df
    
    def create_heatmap(self, output_file):
        """Create heatmap of protein counts"""
        df = self.create_comparison_table()
        
        plt.figure(figsize=(14, 10))
        sns.heatmap(df, annot=True, cmap='YlOrRd', fmt='d',
                   cbar_kws={'label': 'Number of Proteins'})
        plt.title('Secretion System Protein Counts\nComparative Analysis Across Genomes',
                 fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        print(f"✅ Created: {output_file}")
    
    def create_presence_matrix(self, output_file):
        """Create presence/absence matrix"""
        df = self.create_comparison_table()
        presence = df > 0
        
        plt.figure(figsize=(14, 10))
        sns.heatmap(presence, annot=True, cmap='Blues', fmt='d',
                   cbar_kws={'label': 'Present (1) / Absent (0)'})
        plt.title('Secretion System Presence/Absence Matrix\nComparative Analysis Across Genomes',
                 fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        print(f"✅ Created: {output_file}")
    
    def create_summary_report(self, output_file):
        """Create summary report"""
        with open(output_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("SECRETION SYSTEM COMPARATIVE ANALYSIS\n")
            f.write("="*80 + "\n\n")
            
            for genome in self.genomes:
                data = self.data[genome]
                systems = data.get('systems_present', {})
                total = sum(len(p) for p in systems.values())
                
                f.write(f"Genome: {genome}\n")
                f.write(f"  Total proteins: {total}\n")
                f.write(f"  Systems detected: {len([s for s in systems.keys() if systems[s]])}\n")
                f.write(f"  Systems: {', '.join([s for s in systems.keys() if systems[s]])}\n\n")
            
            f.write("="*80 + "\n")
            f.write("COMPARATIVE TABLE\n")
            f.write("="*80 + "\n")
            df = self.create_comparison_table()
            f.write(df.to_string())
        
        print(f"✅ Created: {output_file}")

def main():
    print("="*80)
    print("COMPARATIVE SECRETION SYSTEM ANALYSIS")
    print("="*80)
    
    analyzer = ComparativeAnalyzer()
    analyzer.create_heatmap('comparative_heatmap.png')
    analyzer.create_presence_matrix('presence_matrix.png')
    analyzer.create_summary_report('comparative_summary.txt')

if __name__ == "__main__":
    main()
