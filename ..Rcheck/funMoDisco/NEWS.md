
# funMoDisco 1.0.0

## Initial CRAN Release

### Features

- **Motif Discovery in Functional Data**:
  - **ProbKMA**: Implements a probabilistic K-means algorithm that leverages local alignment and fuzzy clustering to discover recurring patterns (functional motifs) within and across curves.
    - Capable of handling diverse motifs through a family of distances and normalization techniques.
    - Learns motif lengths in a data-driven manner and supports local clustering for misaligned data.
  - **FunBIalign**: Provides hierarchical agglomerative clustering using the Mean Squared Residue Score for motif identification of specified lengths in functional data.
    - Offers a more deterministic approach with user-tunable parameters for control over motif detection.

- **Simulation Tools**: Includes functions to simulate functional data embedded with motifs, enabling users to create benchmark datasets for validating and comparing motif discovery methods.

### Additional Notes
- **Authors**: Marzia Angela Cremona, Francesca Chiaromonte, Jacopo Di Iorio, Niccolo Feresini, Riccardo Lazzarini.
- **License**: GPL (>= 2).
- **System Requirements**: C++20.
