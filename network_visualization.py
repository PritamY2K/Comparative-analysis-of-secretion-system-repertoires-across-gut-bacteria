#!/usr/bin/env python3
"""
Network visualization for secretion systems
"""

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from collections import defaultdict

class NetworkVisualizer:
    def __init__(self):
        self.system_colors = {
            'T1SS': '#FF6B6B', 'T2SS': '#4ECDC4', 'T3SS': '#45B7D1',
            'T4SS': '#96CEB4', 'T5SS': '#FFEAA7', 'T6SS': '#DDA0DD',
            'T8SS': '#FF9FF3', 'T9SS': '#C9B1FF', 'FlgSS': '#54A0FF',
            'T4PSS': '#A8E6CF', 'SecSS': '#FFB347', 'TatSS': '#FD79A8'
        }
        
        self.interactions = {
            'T1SS': [('TolC', 'MacB'), ('TolC', 'LptB'), ('TolC', 'MsbA'),
                    ('MacB', 'BtuC'), ('BtuC', 'BtuD'), ('BtuC', 'BtuF'),
                    ('ZnuA', 'ZnuB'), ('ZnuA', 'ZnuC')],
            'T3SS': [('SsaV', 'SsaT'), ('SsaV', 'SsaS'), ('SsaV', 'SsaR'),
                    ('PrgH', 'PrgI'), ('PrgH', 'PrgJ'), ('PrgH', 'PrgK'),
                    ('InvA', 'InvG')],
            'T6SS': [('Hcp', 'VgrG'), ('Hcp', 'ImpA'), ('VgrG', 'ImpA'),
                    ('ImpA', 'ImpH'), ('VasD', 'VasG')],
            'FlgSS': [('FlgE', 'FlgF'), ('FlgE', 'FlgG'), ('FliC', 'FliD'),
                     ('FliF', 'FliG'), ('FliG', 'FliH'), ('FliM', 'FliN')]
        }
    
    def create_network(self, system, proteins, output_file):
        """Create a single network"""
        G = nx.Graph()
        
        # Add nodes
        for protein in proteins:
            G.add_node(protein)
        
        # Add interactions
        if system in self.interactions:
            for p1, p2 in self.interactions[system]:
                if p1 in proteins and p2 in proteins:
                    G.add_edge(p1, p2)
        
        # Create figure
        plt.figure(figsize=(10, 8))
        
        if len(G.nodes()) > 1:
            pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
        else:
            pos = {proteins[0]: (0, 0)} if proteins else {}
        
        # Draw
        nx.draw_networkx_nodes(G, pos, node_color=self.system_colors.get(system, '#45B7D1'),
                              node_size=2000, alpha=0.8, edgecolors='black', linewidths=2)
        nx.draw_networkx_edges(G, pos, edge_color='gray', width=2, alpha=0.6)
        nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold')
        
        plt.title(f'{system} Protein Interaction Network\n{len(proteins)} proteins, {len(G.edges())} interactions',
                 fontsize=12, fontweight='bold')
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        print(f"✅ Created: {output_file}")

def main():
    # Example usage
    visualizer = NetworkVisualizer()
    
    # Example proteins for each system
    proteins = {
        'T1SS': ['TolC', 'MacB', 'LptB', 'MsbA', 'BtuC', 'BtuD', 'BtuF'],
        'T3SS': ['SsaV', 'SsaT', 'SsaS', 'SsaR', 'PrgH', 'PrgI', 'PrgJ', 'PrgK', 'InvA', 'InvG'],
        'T6SS': ['Hcp', 'VgrG', 'ImpA', 'ImpH', 'VasG', 'VasD'],
        'FlgSS': ['FlgE', 'FlgF', 'FlgG', 'FliC', 'FliD', 'FliF', 'FliG', 'FliH', 'FliM', 'FliN']
    }
    
    for system, prots in proteins.items():
        visualizer.create_network(system, prots, f'network_{system}.png')

if __name__ == "__main__":
    main()
