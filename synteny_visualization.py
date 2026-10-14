#!/usr/bin/env python3
"""
Synteny visualization for secretion systems
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

class SyntenyVisualizer:
    def __init__(self):
        self.system_colors = {
            'T1SS': '#FF6B6B', 'T2SS': '#4ECDC4', 'T3SS': '#45B7D1',
            'T4SS': '#96CEB4', 'T5SS': '#FFEAA7', 'T6SS': '#DDA0DD',
            'T8SS': '#FF9FF3', 'T9SS': '#C9B1FF', 'FlgSS': '#54A0FF',
            'T4PSS': '#A8E6CF', 'SecSS': '#FFB347', 'TatSS': '#FD79A8'
        }
    
    def create_synteny(self, genome_data, output_file):
        """Create synteny plot for a genome"""
        fig, axes = plt.subplots(len(genome_data), 1, figsize=(14, 2 * len(genome_data)))
        fig.suptitle('Synteny of Secretion Systems', fontsize=14, fontweight='bold')
        
        for idx, (system, proteins) in enumerate(genome_data.items()):
            ax = axes[idx] if len(genome_data) > 1 else axes
            
            if not proteins:
                ax.text(0.5, 0.5, f'{system}: Not detected',
                       ha='center', va='center', fontsize=10, color='gray')
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            
            y_pos = 0
            for i, protein in enumerate(proteins[:15]):
                x_pos = i * 0.8
                ax.arrow(x_pos, y_pos, 0.6, 0,
                        head_width=0.12, head_length=0.06,
                        fc=self.system_colors.get(system, '#45B7D1'),
                        ec='black', alpha=0.8, linewidth=1)
                ax.text(x_pos + 0.3, y_pos + 0.15, protein,
                       ha='center', va='bottom', fontsize=7, rotation=45)
            
            ax.text(-0.05, 0.5, system, transform=ax.transAxes,
                   ha='right', va='center', fontsize=10, fontweight='bold')
            
            ax.set_xlim(-0.5, len(proteins[:15]) * 0.8 + 0.5)
            ax.set_ylim(-0.4, 0.5)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_facecolor('#f8f9fa')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        print(f"✅ Created: {output_file
