"""
Molecular SHAP Visualization Module
====================================

Generates publication-quality SHAP importance visualizations for molecules
across multiple model predictions using RDKit and Matplotlib.

Author: Data Science
Date: 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize, LinearSegmentedColormap
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import warnings

warnings.filterwarnings('ignore', category=DeprecationWarning)


class MolecularSHAPVisualizer:
    """
    Generates SHAP importance visualizations for molecules across multiple models.
    
    Attributes:
        highlight_percentile (float): Percentile threshold for consensus highlighting (0-1)
        dpi (int): Resolution for saved figures
        figure_size (Tuple[int, int]): Size of output figures
        font_size (int): Base font size for labels
    """
    
    def __init__(
        self,
        highlight_percentile: float = 0.8,
        dpi: int = 300,
        figure_size: Tuple[int, int] = (10, 8),
        font_size: int = 12
    ):
        """Initialize the visualizer with configuration parameters."""
        self.highlight_percentile = highlight_percentile
        self.dpi = dpi
        self.figure_size = figure_size
        self.font_size = font_size
        self.cmap_reverse = cm.RdYlBu  # Red=high importance, Blue=low
        
    def _mol_from_smiles(self, smiles_string: str) -> Chem.Mol:
        """Parse SMILES and generate 2D coordinates."""
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles_string}")
        
        AllChem.Compute2DCoords(mol)
        return mol
    
    def _get_atom_colors_normalized(
        self,
        shap_vector: np.ndarray
    ) -> Tuple[np.ndarray, Tuple[float, float]]:
        """
        Normalize SHAP values to [0, 1] for colormap and return RGB colors.
        
        Returns:
            colors (np.ndarray): Array of RGB tuples
            value_range (Tuple): (min_value, max_value) for colorbar
        """
        min_val = np.min(shap_vector)
        max_val = np.max(shap_vector)
        
        # Avoid division by zero
        if max_val == min_val:
            normalized = np.ones_like(shap_vector) * 0.5
        else:
            normalized = (shap_vector - min_val) / (max_val - min_val)
        
        # Apply reversed colormap (blue=low, red=high)
        colors = self.cmap_reverse(1 - normalized)  # Reverse for intuitive color
        
        return colors, (min_val, max_val)
    
    def _draw_atom_with_color(
        self,
        mol: Chem.Mol,
        atom_colors: np.ndarray,
        value_range: Tuple[float, float],
        model_idx: int,
        y_pred: float,
        output_path: Path
    ) -> None:
        """Draw molecule with atoms colored by SHAP importance."""
        fig, ax = plt.subplots(figsize=self.figure_size, dpi=self.dpi)
        
        # Get atom positions
        conf = mol.GetConformer()
        pos = {}
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos[idx] = conf.GetAtomPosition(idx)
        
        # Draw using matplotlib
        ax.set_xlim(-2, 12)
        ax.set_ylim(-2, 10)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Draw bonds
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            x = [pos[begin_idx].x, pos[end_idx].x]
            y = [pos[begin_idx].y, pos[end_idx].y]
            ax.plot(x, y, 'k-', linewidth=2, zorder=1)
        
        # Draw atoms with SHAP colors
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            x, y = pos[idx].x, pos[idx].y
            
            # Get color from SHAP value
            color = atom_colors[idx]
            
            # Draw atom circle
            circle = plt.Circle((x, y), 0.35, color=color, ec='black', linewidth=1.5, zorder=2)
            ax.add_patch(circle)
            
            # Draw atom symbol
            ax.text(x, y, symbol, ha='center', va='center', fontsize=self.font_size-2,
                   weight='bold', color='black', zorder=3)
        
        # Add colorbar
        sm = cm.ScalarMappable(cmap=self.cmap_reverse, norm=Normalize(vmin=value_range[0], vmax=value_range[1]))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('SHAP Value', fontsize=self.font_size)
        
        # Title
        title = f"Model {model_idx} | Predicted Y = {y_pred:.3f}"
        ax.set_title(title, fontsize=self.font_size+2, weight='bold', pad=20)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"✓ Saved: {output_path}")
    
    def _get_consensus_highlights(
        self,
        shap_vectors: np.ndarray
    ) -> Tuple[np.ndarray, Dict[int, List[int]]]:
        """
        Identify atoms in top percentile for any model.
        
        Returns:
            highlighted_atoms (np.ndarray): Boolean array of highlighted atoms
            highlight_sources (Dict): {atom_idx: [model_indices]}
        """
        n_models, n_atoms = shap_vectors.shape
        highlighted_atoms = np.zeros(n_atoms, dtype=bool)
        highlight_sources = {i: [] for i in range(n_atoms)}
        
        for model_idx in range(n_models):
            threshold = np.percentile(shap_vectors[model_idx], self.highlight_percentile * 100)
            for atom_idx in range(n_atoms):
                if shap_vectors[model_idx, atom_idx] >= threshold:
                    highlighted_atoms[atom_idx] = True
                    highlight_sources[atom_idx].append(model_idx)
        
        return highlighted_atoms, highlight_sources
    
    def _format_model_labels(self, model_indices: List[int]) -> str:
        """Format model indices as compact string."""
        if len(model_indices) <= 3:
            return ", ".join([f"M{i}" for i in model_indices])
        else:
            # Condense format: M0-2,4
            ranges = []
            start = model_indices[0]
            end = model_indices[0]
            
            for i in model_indices[1:]:
                if i == end + 1:
                    end = i
                else:
                    if start == end:
                        ranges.append(f"M{start}")
                    else:
                        ranges.append(f"M{start}-{end}")
                    start = i
                    end = i
            
            if start == end:
                ranges.append(f"M{start}")
            else:
                ranges.append(f"M{start}-{end}")
            
            return ",".join(ranges)
    
    def _draw_consensus_highlights(
        self,
        mol: Chem.Mol,
        highlighted_atoms: np.ndarray,
        highlight_sources: Dict[int, List[int]],
        output_path: Path
    ) -> None:
        """Draw molecule with consensus highlighted atoms and model labels."""
        fig, ax = plt.subplots(figsize=self.figure_size, dpi=self.dpi)
        
        # Get atom positions
        conf = mol.GetConformer()
        pos = {}
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos[idx] = conf.GetAtomPosition(idx)
        
        # Draw using matplotlib
        ax.set_xlim(-2, 12)
        ax.set_ylim(-2, 10)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Draw bonds
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            x = [pos[begin_idx].x, pos[end_idx].x]
            y = [pos[begin_idx].y, pos[end_idx].y]
            ax.plot(x, y, 'k-', linewidth=2, zorder=1)
        
        # Draw atoms
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            x, y = pos[idx].x, pos[idx].y
            
            # Color based on highlight status
            if highlighted_atoms[idx]:
                color = '#FF6B6B'  # Red for highlighted
                ec_width = 2.5
                ec_color = '#C92A2A'
            else:
                color = '#D0D0D0'  # Light gray for non-highlighted
                ec_width = 1.5
                ec_color = '#808080'
            
            # Draw atom circle
            circle = plt.Circle((x, y), 0.35, color=color, ec=ec_color, 
                               linewidth=ec_width, zorder=2)
            ax.add_patch(circle)
            
            # Draw atom symbol
            ax.text(x, y, symbol, ha='center', va='center', fontsize=self.font_size-2,
                   weight='bold', color='black', zorder=3)
        
        # Draw model labels for highlighted atoms
        for atom_idx, model_indices in highlight_sources.items():
            if highlighted_atoms[atom_idx] and model_indices:
                x, y = pos[atom_idx].x, pos[atom_idx].y
                label = self._format_model_labels(model_indices)
                
                # Offset label above atom
                ax.text(x, y + 0.65, label, ha='center', va='bottom',
                       fontsize=self.font_size-3, weight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', 
                                alpha=0.7, edgecolor='black', linewidth=0.5),
                       zorder=4)
        
        # Title
        ax.set_title("Consensus Highlights - Atoms with Top 20% SHAP in Any Model",
                    fontsize=self.font_size+2, weight='bold', pad=20)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"✓ Saved: {output_path}")
    
    def visualize(
        self,
        y_pred: List[float],
        shap_vectors: np.ndarray,
        smiles_string: str,
        output_dir: str = "./shap_plots"
    ) -> Dict[str, str]:
        """
        Generate SHAP visualizations for molecule across multiple models.
        
        Parameters:
            y_pred (List[float]): Model predictions (length = n_models)
            shap_vectors (np.ndarray): SHAP values, shape (n_models, n_atoms)
            smiles_string (str): SMILES representation of molecule
            output_dir (str): Directory for output images
        
        Returns:
            Dict[str, str]: Paths to generated images
                - "model_{i}" -> path to individual model visualization
                - "consensus" -> path to consensus visualization
        
        Raises:
            ValueError: If SMILES is invalid or array shapes mismatch
        """
        # Validate inputs
        shap_vectors = np.asarray(shap_vectors)
        y_pred = np.asarray(y_pred)
        
        if shap_vectors.ndim != 2:
            raise ValueError(f"shap_vectors must be 2D, got shape {shap_vectors.shape}")
        
        n_models, n_atoms = shap_vectors.shape
        if len(y_pred) != n_models:
            raise ValueError(f"y_pred length ({len(y_pred)}) != n_models ({n_models})")
        
        # Create molecule and output directory
        mol = self._mol_from_smiles(smiles_string)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        output_files = {}
        
        # Generate individual model visualizations
        print(f"\n{'='*60}")
        print(f"Generating SHAP visualizations for {n_models} models")
        print(f"{'='*60}")
        
        for model_idx in range(n_models):
            print(f"\n[Model {model_idx}] Processing...")
            shap_vector = shap_vectors[model_idx]
            atom_colors, value_range = self._get_atom_colors_normalized(shap_vector)
            
            output_file = output_path / f"model_{model_idx}_shap.png"
            self._draw_atom_with_color(
                mol, atom_colors, value_range,
                model_idx, y_pred[model_idx], output_file
            )
            output_files[f"model_{model_idx}"] = str(output_file)
        
        # Generate consensus visualization
        print(f"\n[Consensus] Processing...")
        highlighted_atoms, highlight_sources = self._get_consensus_highlights(shap_vectors)
        consensus_file = output_path / "consensus_highlights.png"
        self._draw_consensus_highlights(mol, highlighted_atoms, highlight_sources, consensus_file)
        output_files["consensus"] = str(consensus_file)
        
        print(f"\n{'='*60}")
        print(f"✓ All visualizations generated successfully!")
        print(f"{'='*60}\n")
        
        return output_files


def visualize_shap_on_molecule(
    y_pred: List[float],
    shap_vectors: np.ndarray,
    smiles_string: str,
    output_dir: str = "./shap_plots",
    highlight_percentile: float = 0.8
) -> Dict[str, str]:
    """
    Generate SHAP visualizations for a molecule across multiple models.
    
    This is the main entry point function. It creates a visualizer instance
    and generates both individual model and consensus visualizations.
    
    Parameters:
        y_pred (List[float]): 1D array of model predictions
            Length = n_models
            Each value represents a model's prediction
        
        shap_vectors (np.ndarray): 2D array of SHAP importance values
            Shape: (n_models, n_atoms_in_molecule)
            shap_vectors[i][j] = SHAP value of atom j for model i
        
        smiles_string (str): SMILES representation of the molecule
            Used to generate 2D molecular structure
        
        output_dir (str, optional): Directory for saving visualizations
            Default: "./shap_plots"
        
        highlight_percentile (float, optional): Percentile for consensus highlighting
            Default: 0.8 (top 20%)
            Valid range: 0.0 - 1.0
    
    Returns:
        Dict[str, str]: Dictionary mapping visualization types to file paths
            Keys:
            - "model_{i}" -> individual model visualization
            - "consensus" -> consensus highlights visualization
    
    Example:
        >>> y_pred = [0.85, 0.92, 0.78]
        >>> shap_vectors = np.array([
        ...     [0.1, 0.8, 0.3, 0.05, 0.9],
        ...     [0.05, 0.7, 0.4, 0.1, 0.85],
        ...     [0.2, 0.6, 0.9, 0.02, 0.7]
        ... ])
        >>> smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
        >>> output = visualize_shap_on_molecule(y_pred, shap_vectors, smiles)
    """
    visualizer = MolecularSHAPVisualizer(highlight_percentile=highlight_percentile)
    return visualizer.visualize(y_pred, shap_vectors, smiles_string, output_dir)
