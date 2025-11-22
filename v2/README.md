# PathWeigh II – Graph-Based Belief Propagation for Pathway Activity Analysis

<img src="./data/image.png" width="400" height="400" alt="Pathway Network Analysis" />

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2025.689915-blue)](https://www.biorxiv.org/)

## Overview

PathWeigh II is an enhanced version of [PathWeigh](https://github.com/zurkin1/PathWeigh), introducing **graph-based network inference** for pathway activity calculation. While PathWeigh uses simple averaging of interaction activities, PathWeigh II employs **NetworkX-based graph decomposition** and **Gaussian-scaled belief propagation** to properly handle pathway topology.

### Key Improvements Over PathWeigh

1. **NetworkX Foundation**: Pathways represented as directed multi-graphs enabling efficient graph algorithms
2. **Component Decomposition**: Weakly connected components analyzed independently for principled handling of disconnected subnetworks
3. **Gaussian Scaling**: Interaction activity penalizes input-output belief inconsistency
4. **Simplified Codebase**: ~40% reduction in code complexity with improved readability and maintainability
5. **Enhanced Performance**: Efficient parallel processing architecture shared with [PathSingle](https://github.com/zurkin1/PathSingle)

### Current Limitations

- **RNA-seq only**: Currently supports RNA-seq data (negative binomial distribution). Microarray support will be added in future releases.

---

## Algorithm

### 1. UDP Calculation

For each gene across samples, fit a probability distribution to expression data:
- **RNA-seq**: Negative binomial distribution
- Calculate UDP: probability of being in "Up" state

### 2. Pathway as Directed Multi-Graph

Transform pathway topology into a NetworkX directed multi-graph:
- **Nodes**: Genes/proteins from pathway
- **Edges**: Interaction relationships (activation/inhibition)
- **Edge attributes**: Store interaction type and parameters

### 3. Gaussian-Scaled Belief Propagation

**Initialization**:
```
beliefs = {gene: UDP[gene] for all genes}
```

**Component Decomposition**:
```python
components = nx.weakly_connected_components(pathway_graph)
```

**Interaction Activity** (for each interaction source → target):
```python
input_belief = max(UDP[source_genes])
output_belief = UDP[target_gene]
activity = input_belief * exp(-(input_belief - output_belief)² / (2σ²))
```

The Gaussian scaling factor:
- Returns 1.0 when input and output beliefs agree
- Decreases exponentially when they disagree
- σ (default: 0.5) controls sensitivity to inconsistency

**Aggregation**:
```
component_activity = mean(interaction activities in component)
pathway_activity = mean(component activities)
```

### Why Gaussian Scaling?

Traditional averaging treats all interactions equally regardless of biological consistency. Gaussian scaling penalizes interactions where upstream activation fails to propagate downstream—capturing post-transcriptional regulation, feedback inhibition, and other biological phenomena that create input-output discordance.

---

## Installation

```bash
git clone https://github.com/zurkin1/PathWeigh.git
cd PathWeigh/v2
pip install -r requirements.txt
```

### Requirements
- Python 3.8+
- pandas
- numpy
- scipy
- scikit-learn
- networkx

---

## Usage

### Basic Workflow

1. **Calculate UDP from RNA-seq data**:
```python
from udp import udp_main

# Prepare input.csv: rows=genes, columns=samples
udp_main()
# Output: data/output_udp.csv
```

2. **Calculate pathway activities**:
```python
from activity import calc_activity

# Uses pathway_relations.csv and output_udp.csv
activity_df = calc_activity()
# Output: data/output_activity.csv
```

### Configuration

Edit constants in `config.py`:

```python
SIGMA = 0.5  # Gaussian scaling parameter (0-1)
```

### Pathway File Format

`pathway_relations.csv`:
```csv
pathway,source,interactiontype,target
p53_pathway,TP53*MDM2,activation,CDKN1A
p53_pathway,ATM,activation$phosphorylation,TP53
p53_pathway,MDM2,inhibition$ubiquitination,TP53
```

- **source/target**: Gene names (use `*` for multiple genes in OR logic)
- **interactiontype**: activation, inhibition, phosphorylation, etc.

---

## Pathway Database

PathWeigh II includes 357 curated pathways from KEGG and BioCarta. View the [full pathway list](https://github.com/zurkin1/PathWeigh/blob/master/v2/data/pathway_relations.csv).

---

## Performance

- **UDP fitting**: ~1-2 minutes for 20,000 genes × 100 samples (parallelized)
- **Activity inference**: ~0.1-1 second per pathway per sample
- **Total runtime**: 357 pathways × 100 samples ≈ 2-5 minutes on 8-core CPU

Parallelization strategies:
- Sample-level: ProcessPoolExecutor across samples
- Gene-level: Multiprocessing for UDP fitting

---

## Benchmark Results

Evaluated on TCGA CRC dataset (577 samples) using MSI status classification:

| Method | Clustering Accuracy |
|--------|---------------------|
| GSEA | 0.443 |
| PathWeigh | 0.433 |
| **PathWeigh II** | **0.635** |

PathWeigh II achieves improved clustering performance through component decomposition and Gaussian scaling, which better captures biologically meaningful pathway perturbations.

---

## Related Projects

- **[PathWeigh](https://github.com/zurkin1/PathWeigh)**: Original pathway activity analysis (Livne & Efroni, 2022)
- **[PathSingle](https://github.com/zurkin1/PathSingle)**: Single-cell pathway analysis

---

## Citation

If you use PathWeigh II in your research, please cite:

```bibtex
@article{livne2025pathweigh2,
  title={PathWeigh II: Graph-Based Belief Propagation for Pathway Activity Analysis},
  author={Livne, Dani},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.689915}
}
```
---

## License

Released under MIT License. See [LICENSE](LICENSE) file.

---

## Contact

**Dani Livne**  
📧 dani.livne@yahoo.com  
GitHub: [@zurkin1](https://github.com/zurkin1)
