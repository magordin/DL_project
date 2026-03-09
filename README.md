    isoform_representation_project/
    │
    ├── README.md
    ├── LICENSE
    ├── environment.yml
    │
    ├── data/
    │   ├── raw/                    # Original datasets (RNA-seq counts, metadata)
    │   │   ├── gene_counts/
    │   │   ├── isoform_counts/
    │   │   └── metadata/
    │   │
    │   ├── processed/              # Cleaned matrices
    │   │   ├── gene_expression/
    │   │   ├── isoform_expression/
    │   │   └── normalized/
    │   │
    │   └── features/               # Feature representations used by models
    │       ├── gene_features/
    │       ├── isoform_features/
    │       ├── exon_features/
    │       └── embeddings/
    │
    ├── notebooks/
    │   ├── exploration_gene_expression.ipynb
    │   ├── exploration_isoform_expression.ipynb
    │   └── representation_comparison.ipynb
    │
    ├── src/
    │   ├── __init__.py
    │   │
    │   ├── config.py               # Paths, experiment parameters
    │   │
    │   ├── data_loading.py         # Load expression matrices
    │   ├── preprocessing.py        # Normalization, filtering
    │   │
    │   ├── representations/        # Core concept of the project
    │   │   ├── gene_representation.py
    │   │   ├── isoform_representation.py
    │   │   ├── exon_representation.py
    │   │   └── embedding_representation.py
    │   │
    │   ├── models/
    │   │   ├── classifier.py
    │   │   ├── neural_network.py
    │   │   └── baseline_models.py
    │   │
    │   ├── evaluation/
    │   │   ├── metrics.py
    │   │   ├── cross_validation.py
    │   │   └── comparison.py       # compare representations
    │   │
    │   ├── plotting.py
    │   └── utils/
    │       ├── io.py
    │       └── helpers.py
    │
    ├── scripts/
    │   ├── 00_submit.sh            # HPC submission
    │   │
    │   ├── 01_preprocess_data.py
    │   ├── 02_build_representations.py
    │   ├── 03_train_models.py
    │   ├── 04_evaluate_models.py
    │   └── 05_plot_results.py
    │
    ├── experiments/                # Each run stored separately
    │   ├── gene_expression/
    │   ├── isoform_expression/
    │   └── combined_features/
    │
    ├── results/
    │   ├── figures/
    │   ├── metrics/
    │   └── logs/
    │
    ├── tests/
    │
    └── docs/
        ├── report/
        ├── slides/
        └── manuscript/