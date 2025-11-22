"""
Benchmark PathBayes on CRC TCGA data comparing against CMSclassifier results.
Classifies colorectal cancer samples into CMS subtypes (CMS1-4) using pathway activities.
Data:
Level 3 TCGA RNA-seq data, downloaded from the TCGA data portal in January 2014.
The data is RSEM-normalized.
After normalization, the data underwent log-transformation (log2 scale).
Non-tumor samples were removed.
Batch effects were corrected using ComBat method.
Additional quality control steps like outlier detection and normalization (including median centering for some data sets) were applied.
"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import Normalizer # Unit norm, row wise. # StandardScaler # Normal distribution. MinMaxScaler # [0,1] range, column wise.
from sklearn.decomposition import PCA
from metrics import *
import decoupler
from pathlib import Path
from scipy.stats import pearsonr
from difflib import get_close_matches


def compare_pathway_results():
        # Load activity outputs
    v1_activity = pd.read_csv('./data/output_activity_v1.csv', index_col=0)
    v2_activity = pd.read_csv('./data/output_activity.csv', index_col=0)

    # Normalize pathway names
    v1_activity.index = v1_activity.index.str.lower().str.strip()
    v2_activity.index = v2_activity.index.str.lower().str.strip()

    # Find common pathways and common samples.
    #common_pathways = list(set(v1_activity.index) & set(v2_activity.index))
    matches_list = []
    for v2_name in list(v2_activity.index):
        matches = get_close_matches(v2_name, list(v1_activity.index), n=1, cutoff=0.6)
        if matches:
            matches_list.append((v2_name, matches[0]))
    common_samples = list(set(v1_activity.columns) & set(v2_activity.columns))

    print(f"Common pathways in activity files: {len(matches_list)}")
    print(f"Common samples: {len(common_samples)}")

    # Subset to common
    v2_names = [match[0] for match in matches_list] # PathWeigh II names
    v1_names = [match[1] for match in matches_list] # PathWeigh v1 names
    v1_common = v1_activity.loc[v1_names, common_samples]
    v2_common = v2_activity.loc[v2_names, common_samples]
    v1_common.index = v2_names # Now both have same index

    # Calculate correlation per sample
    correlations = []
    for sample in tqdm(common_samples):
        r, p = pearsonr(v1_common[sample], v2_common[sample])
        correlations.append(r)

    print(f"\nMean Pearson r across samples: {np.mean(correlations):.3f}")
    print(f"Std: {np.std(correlations):.3f}")
    print(f"Range: {np.min(correlations):.3f} - {np.max(correlations):.3f}")

    # Or: flatten both matrices and compute single correlation
    v1_flat = v1_common.values.flatten()
    v2_flat = v2_common.values.flatten()
    r_overall, p_overall = pearsonr(v1_flat, v2_flat)
    print(f"\nOverall correlation (all pathways × samples): r={r_overall:.3f}, p={p_overall:.2e}")


def load_tcga_data(expression_file, labels_file):
    """Load TCGA CRC expression data and CMS labels"""
    
    # Load activity data (genes × samples)  # (samples (rows) × pathways (columns)) 314 pathways.
    activity = pd.read_csv(expression_file, index_col=0).T
    #print(activity.head())
    #activity = activity.iloc[:100,:]
    print(f"Data shape: {activity.shape}")
    
    # Load CMS labels
    labels = pd.read_csv(labels_file)
    
    # Filter outliers.
    outliers = pd.read_csv('./data/TCGACRC_expression_merged_outlier_strict.txt')
    outliers2 = set(labels.loc[pd.isna(labels["microsatelite"]) | (labels["microsatelite"] == 'Indeterminate')]['id'])
    outset = set(outliers.name).union(outliers2)
    activity = activity.loc[~activity.index.isin(outset)]
    print(f"activity data shape after outlier cleanup: {activity.shape}")
    
    # Filter to TCGA samples only
    #tcga_labels = labels[labels['dataset'] == 'tcga'].copy()
    #tcga_labels = tcga_labels[tcga_labels['CMS_final_network_plus_RFclassifier_in_nonconsensus_samples'] != 'NOLBL']
    
    # Match samples between expression and labels
    common_samples = list(set(activity.index) & set(labels['id']))    
    activity = activity[activity.index.isin(common_samples)]
    labels = labels[labels['id'].isin(common_samples)].copy()
    labels.set_index('id', inplace=True)
    labels = labels.loc[common_samples] # Ensure same order
    labels = labels['microsatelite']
    print(f"activity data shape after TCGA filter: {activity.shape}")

    # Extract labels (final consensus)
    #y = labels_matched['microsatelite'] #CMS_final_network_plus_RFclassifier_in_nonconsensus_samples
    #print(pd.Series(y).value_counts().sort_index())
    
    return activity, labels


def select_cms_relevant_pathways(activity_df):
    """
    Select pathways most relevant to CRC CMS classification.
    Based on biological knowledge of CMS subtypes:
    - CMS1: Immune/MSI, hypermutated
    - CMS2: Canonical/epithelial
    - CMS3: Metabolic
    - CMS4: Mesenchymal/stromal
    """
    activity_df = activity_df.T
    # CMS-relevant pathway keywords
    cms_keywords = [
        # Immune (CMS1)
        'immune', 'inflammation', 'cytokine', 'interferon', 'antigen', 'lymphocyte',
        'T cell', 'B cell', 'NK cell', 'chemokine', 'pd1', 'cytotoxic', 'nkc', 'th1', 'tfh', 'th17', 'treg', 'mdsc', 'complement',
        
        # Canonical/WNT (CMS2)
        'wnt', 'signaling', 'cell cycle', 'proliferation', 'egfr', 'erbb', 'mapk', 'pi3k', 'src', 'jak_stat',
        'caspases', 'proteasome', 'cell_cycle', 'translation', 'notch', 'integrin', 'vegf',
        
        # Metabolic (CMS3)
        'metabolism', 'metabolic', 'glycolysis', 'fatty acid', 'citrate', 'sugar', 'glucose', 'fructose',
        'oxidative phosphorylation', 'amino acid', 'sucrose', 'glactose', 'glutomine', 'glutathione', 'nitrogen',
        'tyrosine', 'fatty_acid', 'arachnoid', 'linoleic',
        
        # Mesenchymal (CMS4)
        'tgf', 'beta', 'adhesion', 'ecm', 'extracellular matrix', 'angiogenesis',
        'focal adhesion', 'integrin',
        
        # General cancer pathways
        'cancer', 'p53', 'apoptosis', 'jak-stat', 'mapk', 'pi3k', 'akt',
        'mtor', 'ras', 'nf-kappa', 'hippo', 'notch', 'hedgehog'
    ]
    
    pathway_names = activity_df.index.str.lower()
    
    # Find pathways matching keywords
    selected_mask = pd.Series([False] * len(pathway_names), index=activity_df.index)
    for keyword in cms_keywords:
        selected_mask |= pathway_names.str.contains(keyword, case=False, na=False)
    
    selected_pathways = activity_df[selected_mask]
    
    print(f"Selected {len(selected_pathways)} CMS-relevant pathways from {len(activity_df)}")
    #print("Sample pathways:")
    #for i, p in enumerate(selected_pathways.index[:10]):
    #    print(f"  {i+1}. {p}")
    
    return selected_pathways.T


labels = ['MSI-H', 'MSS', 'MSI-L'] #'CMS1', 'CMS2', 'CMS3', 'CMS4'


def run_gsea_benchmark():
    #Retrieve gene sets.
    def gmt_to_decoupler(pth: Path) -> pd.DataFrame:
        """
        Parse a gmt file to a decoupler pathway dataframe.
        """
        from itertools import chain, repeat

        pathways = {}

        with Path(pth).open("r") as f:
            for line in f:
                name, _, *genes = line.strip().split("\t")
                pathways[name] = genes

        return pd.DataFrame.from_records(
            chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
            columns=["geneset", "genesymbol"],
        )

    reactome = gmt_to_decoupler("./data/c2.cp.reactome.v7.5.1.symbols.gmt")
    reactome = reactome.rename(columns={'geneset': 'source', 'genesymbol': 'target'})
    print(reactome.columns)
    print(reactome.head())
 
    # Load data
    df = pd.read_csv('./data/TCGACRC_expression-merged.zip', sep='\t', index_col=0).T #, nrows=100, usecols=range(100)
    print(df.head())
    labels = pd.read_csv('./data/TCGACRC_clinical-merged.csv')
    # Match samples between expression and labels
    common_samples = list(set(df.index) & set(labels['id']))    
    df = df[df.index.isin(common_samples)]
    labels = labels[labels['id'].isin(common_samples)].copy()
    labels.set_index('id', inplace=True)
    labels = labels.loc[common_samples] # Ensure same order
    labels = labels['microsatelite']
    print(f"data shape after TCGA filter: {df.shape}")

    # decoupler expects the genes as columns.
    # Run GSEA once on all samples
    activity, pvals = decoupler.mt.gsea(
        data=df,
        net=reactome,
        tmin=5,
        verbose=True
    )
    '''
    #On older version of decoupler.
    activity, norm, pvals = decoupler.run_gsea(
    df,
    reactome,
    verbose=True
    )
    '''
    #print(activity.head())
    # Cluster requires (n_samples, n_features).
    # Convert labels to numeric for metrics
    kmeans = clustering(activity, n_clusters=3)
    y_pred_kmeans = kmeans.labels_
    
    # Metrics
    label_map = {label: index for index, label in enumerate(labels)}
    #index_map = {index: label for index, label in enumerate(labels)}
    y_true_numeric = np.array([label_map[label] for label in labels.values])
    silhouette, calinski, special_acc, completeness, homogeneity, adjusted_mi = calc_stats(
        activity, 
        y_true_numeric, 
        y_pred_kmeans,
        debug=True
    )
    activity.to_csv('activity_gsea.csv')


if __name__ == "__main__":
    """Main benchmark pipeline"""
    if True:
        #run_gsea_benchmark()
        compare_pathway_results()
        exit(0)

    # Load data
    activity, y_true = load_tcga_data('./data/output_activity.csv', labels_file='./data/TCGACRC_clinical-merged.csv')

    # Select CMS-relevant pathways (optional - comment out to use all pathways)
    #activity = select_cms_relevant_pathways(activity)
    
    # Scale the data.
    scaler = Normalizer()
    activity = scaler.fit_transform(activity)
    
    # Dimensionality reduction with PCA
    #pca = PCA(n_components=30, svd_solver='arpack')
    #activity = pca.fit_transform(activity)
    #print(f"PCA explained variance ratio: {pca.explained_variance_ratio_.sum():.3f}")
       
    # Cluster requires (n_samples, n_features).
    #print("Clustering (k=3)...")
    # Convert labels to numeric for metrics
    kmeans = clustering(activity, n_clusters=3)
    y_pred_kmeans = kmeans.labels_
    
    # Metrics
    label_map = {label: index for index, label in enumerate(labels)}
    #index_map = {index: label for index, label in enumerate(labels)}
    y_true_numeric = np.array([label_map[label] for label in y_true.values])
    silhouette, calinski, special_acc, completeness, homogeneity, adjusted_mi = calc_stats(
        activity, 
        y_true_numeric, 
        y_pred_kmeans,
        debug=True
    )