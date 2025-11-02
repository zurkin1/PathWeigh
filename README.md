# PathWeigh - Quantifying the Behavior of Pathway Cascades

<figure>
  <img src="https://norbis.w.uib.no/files/2016/05/F1.large_-768x623.jpg" alt="Pathway analysis" />
  <figcaption><small>(Source: https://norbis.w.uib.no/files/2016/05/F1.large_-768x623.jpg)</small></figcaption>
</figure>

## Overview
PathWeigh is an advanced Python-based pathway analysis tool specifically designed for microarray and RNA-sequencing data analysis. This innovative software implements a novel graph-based algorithm, enabling comprehensive analysis of diverse cellular states, including but not limited to T cell subtypes. 

Key attributes of PathWeigh include:

- Open-source architecture
- Extensible framework
- Computational efficiency

These features collectively provide researchers with a robust and versatile platform for extracting biologically meaningful insights from high-dimensional transcriptomics data.

## Key Features
- Tailored for microarray and RNA-seq data analysis
- Unique graph-based algorithm for pathway analysis
- Efficient classification of diverse cellular states
- Open-source and extensible
- Computationally efficient

## Installation
Simply clone this repository using git clone command.

## Input File Requirements
PathWeigh requires a single input file.

- `input.csv`: A gene expression matrix representing cell populations, with dimensions `(n_genes) x (k_cells)`.

Data Specifications:

1. Gene Symbols: 
   - Must be provided in the first column.
   - Can be in either lower or upper case.

2. Expression Data:
   - Should be in non-log scale.
   - An anti-log transformation is automatically applied if the maximum expression value is < 20.

3. Gene Selection:
   - PathWeigh employs internal marker gene selection algorithms.
   - Consequently, not all provided signature genes may be utilized in the analysis.

4. RNA-seq Normalization:
   - Accepted normalization methods: TPM, RPKM, or FPKM.
   - This is crucial as PathWeigh's default behavior assumes a negative binomial distribution of gene expression.

5. Gene Nomenclature:
   - Gene names should adhere to the HUGO standard (e.g., as defined in https://www.genenames.org/).
   - Generally, gene names used should correspond to those in the pathway definitions.

Adherence to these specifications ensures proper formatting of input data, facilitating accurate analysis by PathWeigh.

## Quick Start
```
df = pd.read_csv('./data/input.csv')
calc_activity(df)
```

## Usage
For detailed usage instructions, please refer to the notebooks in the `code` folder, which demonstrate RNA-seq data processing using PathWeigh.

## Supported Pathways
PathWeigh currently supports 581 curated pathways. Click the link to view the full list. [List of supported pathways.](code/data/pathnames.txt)

[Guide for adding a new pathway.](code/data/guide.md)

## Contributing
We welcome contributions! Please see our Contributing Guidelines for more information on how to get involved.

## License
PathWeigh is available under the MIT license. See the LICENSE file for more details.

## Support
For questions, issues, or feature requests, please open an issue on our GitHub repository.

For additional support, contact: zurkin at yahoo dot com.

## Citation
If you use PathSingle in your research, please cite our paper:

Livne, D., Efroni, S. (2022). PathWeigh – Quantifying the Behavior of Biochemical Pathway Cascades. In: Rojas, I., Valenzuela, O., Rojas, F., Herrera, L.J., Ortuño, F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2022. Lecture Notes in Computer Science(), vol 13347. Springer, Cham. https://doi.org/10.1007/978-3-031-07802-6_29

## Acknowledgments
We thank the reviewers of this work for their valuable feedback and contributions.
