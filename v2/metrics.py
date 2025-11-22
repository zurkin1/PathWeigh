import numpy as np
from sklearn.metrics import silhouette_score, calinski_harabasz_score, completeness_score, homogeneity_score, adjusted_mutual_info_score, rand_score, normalized_mutual_info_score
from sklearn.metrics import accuracy_score, confusion_matrix, cohen_kappa_score
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import linear_sum_assignment
from tqdm.notebook import tqdm
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import VarianceThreshold
import scipy
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from umap import UMAP
from tabulate import tabulate
import sys
import os


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(CURRENT_DIR)
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)


def to_prob_logscale(x, vmax=22.0, eps=1e-6):
    p = np.log1p(x) / np.log1p(vmax) # log(1+x)/log(1+vmax)
    return np.clip(p, eps, 1 - eps)


def to_prob_power(x, vmax=22.0, alpha=0.5, eps=1e-6):
    p = (np.power(x, alpha)) / (np.power(vmax, alpha))
    return np.clip(p, eps, 1 - eps)

#General metric functions used in the benchmarks.
#Dunn index calculation function.
def dunn_index(X, labels):
    #Compute pairwise distances.
    distances = squareform(pdist(X))
    
    #Find unique cluster labels.
    unique_labels = np.unique(labels)
    
    #Calculate inter-cluster distances (minimum distance between points in different clusters).
    inter_cluster_dists = np.inf
    for i in range(len(unique_labels)):
        for j in range(i + 1, len(unique_labels)):
            #cluster_i = X[labels == unique_labels[i]]
            #cluster_j = X[labels == unique_labels[j]]
            inter_cluster_dists = min(inter_cluster_dists, np.min(distances[np.ix_(labels == unique_labels[i], labels == unique_labels[j])]))
    
    #Calculate intra-cluster distances (maximum distance within a cluster).
    intra_cluster_dists = 0
    for label in unique_labels:
        cluster = X[labels == label]
        intra_cluster_dists = max(intra_cluster_dists, np.max(pdist(cluster)) if len(cluster) > 1 else 0)
    
    #Calculate Dunn index.
    return inter_cluster_dists / intra_cluster_dists if intra_cluster_dists != 0 else np.inf


#Special accuracy function for clustering.
def acc(y_true, y_pred):
    """
    Calculate clustering accuracy. This acc function calculates clustering accuracy, which is a measure of how well a clustering algorithm assigns
    data points to clusters compared to the ground truth labels. It does so by finding the best matching between the predicted cluster assignments
    (y_pred) and the true labels (y_true), and then calculating the accuracy based on this matching. Here's how the function works:
    1. It first ensures that both y_true and y_pred are integer arrays.
    2. It asserts that the sizes of y_true and y_pred are equal.
    3. It initializes a confusion matrix w where each entry w[i, j] denotes the number of data points that belong to cluster i according to y_pred
       and to cluster j according to y_true.
    4. It computes the optimal one-to-one matching between the predicted clusters and the true clusters using the Hungarian algorithm, implemented
       in linear_assignment from scikit-learn.
    5. It calculates the accuracy by summing the values of the confusion matrix corresponding to the optimal matching and dividing by the total 
       number of data points.
    6. This accuracy measure is specifically tailored for clustering algorithms, where there is no inherent ordering of cluster labels. It 
       accounts for the fact that clusters may be assigned arbitrary labels by the algorithm, and it finds the best matching between these labels and the true labels to compute accuracy.

    Traditional accuracy metrics assume a one-to-one correspondence between predicted and true labels, which is not the case in clustering. 
    Instead, the acc function finds the best matching between clusters and true labels, providing a more appropriate measure of clustering 
    performance. We need to find the best matching between the predicted clusters and the true clusters and then evaluate the performance based on
    this matching. This is what the acc function does by finding the optimal matching using the Hungarian algorithm and computing the accuracy
    based on this matching.
    Arguments
        y: true labels, numpy.array with shape `(n_samples,)`
        y_pred: predicted labels, numpy.array with shape `(n_samples,)`
    Return
        accuracy, in [0,1]
    """
    # Ensure y_true and y_pred are NumPy arrays of integers.
    y_true = np.asarray(y_true, dtype=np.int64)
    y_pred = np.asarray(y_pred, dtype=np.int64)

    assert y_pred.size == y_true.size

    # Determine the size of the confusion matrix and initialize it.
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)

    # Fill the confusion matrix.
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    
    # Perform the assignment using the Hungarian algorithm.
    ind = linear_sum_assignment(w.max() - w)

    # Return accuracy.
    return sum(w[ind[0], ind[1]]) * 1.0 / y_pred.size


def calc_stats(act_mat, true_labels, pred_labels, debug=False):
    # Overall accuracy
    #print(f'overall_acc: {accuracy_score(true_labels, pred_labels)}')
    
    # Confusion matrix
    #cm = confusion_matrix(true_labels, pred_labels)
    #print(f"confusion matrix")
    #df_cm = pd.DataFrame(cm, index=labels, columns=labels)
    #print(cm)

    # Per-class sensitivity (recall)
    #class_acc = {}
    #for i, cms in enumerate(labels):
    #    if cm[i, :].sum() > 0:
    #        class_acc[cms] = cm[i, i] / cm[i, :].sum()
    #    else:
    #        class_acc[cms] = 0.0
    #    print(f'class: {cms} acc: {class_acc[cms]}')
    
    # Cohen's Kappa
    #print(f'Cohen kappa: {cohen_kappa_score(true_labels, pred_labels):.2f}')

    #Silhouette score.
    Silhouette = silhouette_score(act_mat, pred_labels, metric='euclidean')

    #Calinski-Harabasz index.
    Calinski = calinski_harabasz_score(act_mat, pred_labels)

    #Special accuracy function.
    Special = acc(true_labels, pred_labels)

    #Completeness score.
    Completeness =  completeness_score(true_labels, pred_labels)

    #Homogeneity_score.
    Homogeneity = homogeneity_score(true_labels, pred_labels)

    #Adjusted_mutual_info_scor
    Adjusted = adjusted_mutual_info_score(true_labels, pred_labels)
    rs = rand_score(true_labels, pred_labels)
    nmi = normalized_mutual_info_score(true_labels, pred_labels)

    if debug:
        data = [
            ("Silhouette Score", Silhouette),
            ("Calinski-Harabasz Index", Calinski),
            ("Special accuracy", Special),
            ("Completeness score", Completeness),
            ("Homogeneity score", Homogeneity),
            ("Adjusted mutual info score", Adjusted),
            ("Rand score", rs),
            ("Normalized mutual info score", nmi),
        ]
        print(tabulate(data, headers=["Metric", "Value"], tablefmt="github", floatfmt=".6f"))
    '''
    # Plot: True vs Predicted clusters (scatter plot using PCA for visualization)
    embedding = UMAP(n_neighbors=15, min_dist=0.3, random_state=42, n_jobs=1).fit_transform(act_mat)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for i, (labels, title) in enumerate([(true_labels, "True Labels"),
                                        (pred_labels, "Predicted Clusters")]):
        ax[i].scatter(embedding[:, 0], embedding[:, 1],
                    c=labels, cmap='tab10', s=40, alpha=0.7)
        ax[i].set_title(title)
        ax[i].set_aspect('equal', 'box')
    plt.tight_layout()
    plt.show()
    '''
    return Silhouette, Calinski, Special, Completeness, Homogeneity, Adjusted


def feature_importance(scores_df, true_labels):
        # Loop through the features, omitting one at a time and calculating the acc score.
        # Dictionary to store the accuracy results.
        acc_results = {}

        # Loop through each gene set.
        for pathway in tqdm(scores_df.columns):
            # Omit the current pathway.
            pathsingle_results_matrix_omitted = scores_df.drop(columns=[pathway]).values

            selector = VarianceThreshold(threshold=0.01)
            pathsingle_results_matrix_omitted = selector.fit_transform(pathsingle_results_matrix_omitted)

            #Scale the data.
            scaler = MinMaxScaler()
            pathsingle_results_matrix_omitted = scaler.fit_transform(pathsingle_results_matrix_omitted)
            
            # Perform KMeans clustering.
            num_clusters = len(np.unique(true_labels))
            kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(pathsingle_results_matrix_omitted)
            cluster_assignments = kmeans.labels_
            
            # Evaluate clustering results using acc metric.
            clustering_accuracy = acc(true_labels, cluster_assignments)
            
            # Store the accuracy result.
            acc_results[pathway] = clustering_accuracy

        # Sort the gene sets according to the acc result (the best gene set is the one that reduced the results the most).
        sorted_genesets = sorted(acc_results.items(), key=lambda x: x[1])

        # Display the sorted gene sets.
        print("Top 40 gene sets by clustering accuracy:")
        for pathway, accuracy in sorted_genesets[:40]:
            print(f'Pathway set: {pathway}')

def clustering(results_matrix, n_clusters=3):
    #Perform KMeans clustering. (n_samples, n_features)
    model = KMeans(n_clusters, random_state=42).fit(results_matrix)
    #model = SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors', random_state=42).fit(results_matrix)
    return model

def gaussian_scaling(p, q, sigma=0.5):
    """Calculate the scaled activity of inputs (p) of an interaction by the outputs (q).
    1. Distance term: (p - q)**2 : Measures squared difference between p and q. Always positive. Larger when p and q are far apart.
    2. Scaling factor: exp(-(p-q)²/(2σ²)) : Returns 1.0 when p=q. Decreases exponentially as |p-q| increases. σ controls how quickly scaling drops off.
    3. Final value: p * scaling : When p=q: returns p (scaling=1). When p≠q: reduces p based on distance. Never increases above p.
    """
    scaling = np.exp(-(p - q)**2 / (2*sigma**2))
    return p * scaling

def consistency_scaling(p, q):
    """
    Calculate the scaled activity of inputs (p) of an interaction by the outputs (q).
    measures similarity or consistency between p and q, producing values from −1 to 1,
    positive when the pair agrees (both high or both low), negative when they disagree."""
    return p * q - (1 - p) * (1 - q)

def affine_consistency_scaling(p, q):
    '''Scale to [0,1].'''
    return (consistency_scaling(p, q) + 1)/2

def scaled_difference(u, v):
    '''Scaled difference measure highlighting disagreement.'''
    return u*(1-v) + (1-u)*v

def proximity_scaling(p, q):
    """Calculate the scaled activity of inputs (p) of an interaction by the outputs (q)."""
    proximity = 1 - abs(p -q) #How close p and q are.
    return p * proximity

def calculate_sparsity(adata):
    """Calculate data sparsity percentage."""
    if scipy.sparse.issparse(adata.X):
        sparsity2 = (1.0 - (adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))) * 100
    else:
        sparsity2 = np.sum(adata.X == 0) / adata.X.size * 100

    print(f'nnz-based sparsity: {sparsity2:.2f}%.')    
    return sparsity2

def choose_scaling_method(sparsity):
    """Choose scaling method based on data sparsity.
    adata.X.nnz gives number of non-zero elements.
    adata.X.shape[0] * adata.X.shape[1] gives total matrix size.
    """
    return proximity_scaling if sparsity >= 40 else gaussian_scaling