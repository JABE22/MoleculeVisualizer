#!/usr/bin/env python3
"""
Molecular Visualizer - Beautiful Command Line Tool
Convert SMILES strings to publication-quality molecular structure images.

Features:
  - Single molecule visualization
  - Batch processing
  - Stdin support
  - CPK element coloring (standard chemistry colors)
  - Publication-quality output
  - Fixed image dimensions (600x400px default)
  - Professional typography (10pt font)
  - Anti-aliased rendering
  - Beautiful bond visualization

Installation:
  pip install rdkit matplotlib

Quick Start:
  python mol_visualizer.py "CC(=O)OC1=CC=CC=C1C(=O)O"
  python mol_visualizer.py "SMILES" -o output.png
  python mol_visualizer.py --batch molecules.txt -d ./images

Usage:
  Single molecule:    python mol_visualizer.py "SMILES" [options]
  Batch file:        python mol_visualizer.py --batch FILE -d DIR
  Stdin:             cat smiles.txt | python mol_visualizer.py --stdin
"""

import sys
import argparse
from pathlib import Path
from typing import Optional, Tuple, Dict, List
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Circle, Wedge
from matplotlib.font_manager import FontProperties
from rdkit import Chem
from rdkit.Chem import AllChem

warnings.filterwarnings('ignore', category=DeprecationWarning)


class MoleculeVisualizer:
    """
    Beautiful molecular structure visualizer with professional typography.
    
    Features:
      - Fixed image dimensions (600x400 pixels default)
      - CPK color scheme (standard chemistry colors)
      - Professional 10pt typography
      - Anti-aliased rendering
      - Smart bond visualization (single, double, triple)
      - Molecular formula auto-calculation
      - Publication-ready output (300 DPI)
    """
    
    # CPK Color Scheme (standard chemistry colors)
    CPK_COLORS = {
        'H': '#FFFFFF',      # White
        'C': '#909090',      # Gray
        'N': '#3050F8',      # Blue
        'O': '#FF0D0D',      # Red
        'F': '#90E050',      # Light Green
        'P': '#FF8000',      # Orange
        'S': '#FFFF30',      # Yellow
        'Cl': '#1FF01F',     # Green
        'Br': '#A62929',     # Brown
        'I': '#940094',      # Purple
        'Default': '#C0C0C0' # Silver
    }
    
    def __init__(
        self,
        width: int = 600,
        height: int = 400,
        dpi: int = 1200,
        font_size: int = 10
    ):
        """
        Initialize molecular visualizer.
        
        Args:
            width: Image width in pixels (default: 600)
            height: Image height in pixels (default: 400)
            dpi: Resolution in DPI (default: 300, publication quality)
            font_size: Font size in points (default: 10)
        """
        self.width = width
        self.height = height
        self.dpi = dpi
        self.font_size = font_size
        
        # Convert pixels to inches (for matplotlib)
        self.fig_width_inches = width / dpi
        self.fig_height_inches = height / dpi
        
        # Font configuration
        self.title_font = FontProperties(
            family='sans-serif',
            size=font_size,
            weight='normal'
        )
        self.formula_font = FontProperties(
            family='monospace',
            size=font_size - 2,
            weight='normal'
        )
    
    def _smiles_to_mol(self, smiles: str) -> Optional[Chem.Mol]:
        """
        Parse SMILES string and compute 2D coordinates.
        
        Args:
            smiles: SMILES string (e.g., "CC(=O)OC1=CC=CC=C1C(=O)O")
            
        Returns:
            RDKit molecule object or None if invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Add hydrogens for proper visualization
        mol = Chem.AddHs(mol)
        
        # Compute 2D layout
        AllChem.Compute2DCoords(mol)
        
        return mol
    
    def _get_atom_color(self, atom: Chem.Atom) -> str:
        """
        Get CPK color for atom based on element.
        
        Args:
            atom: RDKit atom object
            
        Returns:
            Hex color code (e.g., '#FF0D0D')
        """
        symbol = atom.GetSymbol()
        return self.CPK_COLORS.get(symbol, self.CPK_COLORS['Default'])
    
    def _get_molecular_formula(self, mol: Chem.Mol) -> str:
        """
        Calculate molecular formula from molecule.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Molecular formula string (e.g., "C2H6O")
        """
        from collections import defaultdict
        
        formula_dict = defaultdict(int)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            formula_dict[symbol] += 1
        
        # Sort by common order: C, H, then alphabetical
        priority = ['C', 'H']
        formula = ''
        
        for element in priority:
            if element in formula_dict:
                count = formula_dict[element]
                formula += element if count == 1 else f"{element}{count}"
                del formula_dict[element]
        
        for element in sorted(formula_dict.keys()):
            count = formula_dict[element]
            formula += element if count == 1 else f"{element}{count}"
        
        return formula
    
    def _draw_molecule(self, mol: Chem.Mol, ax):
        """
        Draw molecule structure with bonds and atoms.
        
        Args:
            mol: RDKit molecule object
            ax: Matplotlib axes object
        """
        # Get conformer (2D coordinates)
        conf = mol.GetConformer()
        
        # Draw bonds first (so they appear behind atoms)
        for bond in mol.GetBonds():
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            
            # Get atom positions
            pos_begin = conf.GetAtomPosition(begin_atom_idx)
            pos_end = conf.GetAtomPosition(end_atom_idx)
            
            x = [pos_begin.x, pos_end.x]
            y = [pos_begin.y, pos_end.y]
            
            # Determine bond type and width
            bond_type = bond.GetBondType()
            if bond_type == Chem.BondType.SINGLE:
                linewidth = 1.2
                ax.plot(x, y, 'k-', linewidth=linewidth, solid_capstyle='round')
            elif bond_type == Chem.BondType.DOUBLE:
                # Double bond: two parallel lines
                linewidth = 1.0
                # Perpendicular offset
                dx = pos_end.x - pos_begin.x
                dy = pos_end.y - pos_begin.y
                length = np.sqrt(dx**2 + dy**2)
                
                if length > 0:
                    # Unit perpendicular vector
                    px = -dy / length * 0.08
                    py = dx / length * 0.08
                    
                    # First line
                    ax.plot(
                        [x[0] + px, x[1] + px],
                        [y[0] + py, y[1] + py],
                        'k-', linewidth=linewidth, solid_capstyle='round'
                    )
                    # Second line
                    ax.plot(
                        [x[0] - px, x[1] - px],
                        [y[0] - py, y[1] - py],
                        'k-', linewidth=linewidth, solid_capstyle='round'
                    )
            elif bond_type == Chem.BondType.TRIPLE:
                # Triple bond: three lines
                linewidth = 0.8
                dx = pos_end.x - pos_begin.x
                dy = pos_end.y - pos_begin.y
                length = np.sqrt(dx**2 + dy**2)
                
                if length > 0:
                    px = -dy / length * 0.12
                    py = dx / length * 0.12
                    
                    # Center line
                    ax.plot(x, y, 'k-', linewidth=linewidth, solid_capstyle='round')
                    # Side lines
                    ax.plot(
                        [x[0] + px, x[1] + px],
                        [y[0] + py, y[1] + py],
                        'k-', linewidth=linewidth, solid_capstyle='round'
                    )
                    ax.plot(
                        [x[0] - px, x[1] - px],
                        [y[0] - py, y[1] - py],
                        'k-', linewidth=linewidth, solid_capstyle='round'
                    )
            else:
                # Aromatic or other
                ax.plot(x, y, 'k-', linewidth=1.2, solid_capstyle='round')
        
        # Draw atoms
        for atom in mol.GetAtoms():
            atom_idx = atom.GetIdx()
            pos = conf.GetAtomPosition(atom_idx)
            
            # Skip hydrogen atoms (already implicit in structure)
            if atom.GetSymbol() == 'H':
                continue
            
            # Atom color
            color = self._get_atom_color(atom)
            
            # Draw atom circle with border
            circle = Circle(
                (pos.x, pos.y),
                radius=0.35,
                color=color,
                ec='black',
                linewidth=1,
                zorder=10,
                alpha=0.95
            )
            ax.add_patch(circle)
            
            # Draw element symbol on top of the atom
            label_offset_y = 0.0  # Distance above atom center
            ax.text(
                pos.x,
                pos.y + label_offset_y,
                atom.GetSymbol(),
                ha='center',
                va='center',
                fontsize=6,  # Smaller font for labels
                fontweight='normal',
                color='black',
                zorder=11,
                family='sans-serif'
            )
    
    def draw(
        self,
        smiles: str,
        output_path: str,
        show_formula: bool = True
    ) -> str:
        """
        Render molecule and save as PNG.
        
        Args:
            smiles: SMILES string
            output_path: Output PNG file path
            show_formula: Whether to display molecular formula
            
        Returns:
            Output file path
            
        Raises:
            ValueError: If SMILES is invalid
        """
        # Parse SMILES
        mol = self._smiles_to_mol(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Get molecular formula
        formula = self._get_molecular_formula(mol)
        
        # Create figure with fixed dimensions
        fig, ax = plt.subplots(
            figsize=(self.fig_width_inches, self.fig_height_inches),
            dpi=self.dpi
        )
        
        # Configure axes
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Draw molecule
        self._draw_molecule(mol, ax)
        
        # Set axis limits with padding
        conf = mol.GetConformer()
        xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
        ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
        
        if xs and ys:
            x_min, x_max = min(xs), max(xs)
            y_min, y_max = min(ys), max(ys)
            
            # Add padding
            padding = 1.2
            x_range = x_max - x_min or 1
            y_range = y_max - y_min or 1
            
            ax.set_xlim(x_min - padding, x_max + padding)
            ax.set_ylim(y_min - padding, y_max + padding)
        
        # Add title with SMILES and formula
        title_text = f"{smiles}"
        if show_formula:
            title_text += f"  |  {formula}"
        
        ax.text(
            0.5,
            1.1,
            title_text,
            transform=ax.transAxes,
            ha='center',
            va='top',
            fontsize=5,
            fontproperties=self.title_font,
            bbox=dict(
                facecolor='#F8F8F8',
                edgecolor='#E0E0E0',
                linewidth=0.5,
                alpha=0.85
            )
        )
        
        # Add subtle background
        fig.patch.set_facecolor('white')
        ax.patch.set_facecolor('white')
        
        # Tight layout
        plt.subplots_adjust(
            left=0.05,
            right=0.95,
            top=0.90,
            bottom=0.05
        )
        
        # Save with high quality
        plt.savefig(
            output_path,
            dpi=self.dpi,
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none',
            format='png'
        )
        plt.close(fig)
        
        return output_path


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Convert SMILES strings to beautiful molecule images',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Single molecule:
    python mol_visualizer.py "CC(=O)OC1=CC=CC=C1C(=O)O"
    python mol_visualizer.py "SMILES" -o aspirin.png
  
  Batch processing:
    python mol_visualizer.py --batch molecules.txt -d ./images
  
  Stdin:
    cat smiles.txt | python mol_visualizer.py --stdin
        """
    )
    
    # Input mode
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument(
        'smiles',
        nargs='?',
        help='SMILES string to visualize'
    )
    input_group.add_argument(
        '--batch', '-b',
        metavar='FILE',
        help='Batch file (one SMILES per line)'
    )
    input_group.add_argument(
        '--stdin',
        action='store_true',
        help='Read SMILES from stdin'
    )
    
    # Output options
    parser.add_argument(
        '-o', '--output',
        metavar='FILE',
        help='Output PNG filename (single mode)'
    )
    parser.add_argument(
        '-d', '--output-dir',
        metavar='DIR',
        default='.',
        help='Output directory (default: current directory)'
    )
    
    # Visualization options
    parser.add_argument(
        '--width',
        type=int,
        default=600,
        help='Image width in pixels (default: 600)'
    )
    parser.add_argument(
        '--height',
        type=int,
        default=400,
        help='Image height in pixels (default: 400)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Resolution in DPI (default: 300)'
    )
    parser.add_argument(
        '--font-size',
        type=int,
        default=10,
        help='Font size in points (default: 10)'
    )
    parser.add_argument(
        '--no-formula',
        action='store_true',
        help='Hide molecular formula'
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress output messages'
    )
    
    args = parser.parse_args()
    
    # Determine input mode
    if not args.smiles and not args.batch and not args.stdin:
        parser.print_help()
        sys.exit(1)
    
    # Initialize visualizer
    viz = MoleculeVisualizer(
        width=args.width,
        height=args.height,
        dpi=args.dpi,
        font_size=args.font_size
    )
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process based on mode
    if args.smiles:
        # Single SMILES mode
        try:
            output_filename = args.output or f"molecule_{args.smiles[:30]}.png"
            output_path = output_dir / output_filename
            
            viz.draw(
                args.smiles,
                str(output_path),
                show_formula=not args.no_formula
            )
            
            if not args.quiet:
                print(f"✓ Saved: {output_path}")
            
            sys.exit(0)
        except ValueError as e:
            print(f"✗ Error: {e}", file=sys.stderr)
            sys.exit(1)
    
    elif args.batch:
        # Batch file mode
        batch_file = Path(args.batch)
        if not batch_file.exists():
            print(f"✗ File not found: {batch_file}", file=sys.stderr)
            sys.exit(1)
        
        success_count = 0
        fail_count = 0
        
        with open(batch_file, 'r') as f:
            for idx, line in enumerate(f, 1):
                # Skip empty lines and comments
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    output_filename = f"mol_{idx:04d}.png"
                    output_path = output_dir / output_filename
                    
                    viz.draw(
                        line,
                        str(output_path),
                        show_formula=not args.no_formula
                    )
                    
                    if not args.quiet:
                        print(f"✓ [{idx}] {line[:50]}")
                    
                    success_count += 1
                except ValueError as e:
                    fail_count += 1
                    if not args.quiet:
                        print(f"✗ [{idx}] Invalid SMILES: {line[:50]}")
        
        if not args.quiet:
            print(f"\nResults: {success_count} success, {fail_count} failed")
        
        sys.exit(0 if fail_count == 0 else 1)
    
    elif args.stdin:
        # Stdin mode
        success_count = 0
        fail_count = 0
        
        for idx, line in enumerate(sys.stdin, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            try:
                output_filename = f"mol_{idx:04d}.png"
                output_path = output_dir / output_filename
                
                viz.draw(
                    line,
                    str(output_path),
                    show_formula=not args.no_formula
                )
                
                if not args.quiet:
                    print(f"✓ [{idx}] {line[:50]}")
                
                success_count += 1
            except ValueError:
                fail_count += 1
                if not args.quiet:
                    print(f"✗ [{idx}] Invalid SMILES: {line[:50]}")
        
        if not args.quiet:
            print(f"\nResults: {success_count} success, {fail_count} failed")
        
        sys.exit(0 if fail_count == 0 else 1)


if __name__ == '__main__':
    main()
