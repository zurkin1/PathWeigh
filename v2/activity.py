import pandas as pd
import numpy as np
from metrics import *
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
from graphs import *


def is_inhibitory(interaction_type):
    """Check if interaction type is inhibitory"""
    inhibitory_keywords = ['inhibition', 'repression', 'dissociation', 'dephosphorylation', 'ubiquitination']
    return any(keyword in interaction_type.lower() for keyword in inhibitory_keywords)


def pathway_activity(G, sample_udp):
    """
    Pathway activity using node belief: gaussian scaling of sum(incoming UDP) and sum(outgoing UDP)
    Returns average over nodes.
    """
    node_beliefs = {}
    for node in G.nodes:
        src_genes = G.nodes[node]['source_genes']
        tgt_genes = G.nodes[node]['target_genes']
        itype = G.nodes[node]['interaction_type']
        
        # Compute source and target activities from sample
        src_activity = sum(sample_udp.get(g, 0.0) for g in src_genes)
        tgt_activity = sum(sample_udp.get(g, 0.0) for g in tgt_genes)
        tgt_activity = max(1e-10, tgt_activity)
        
        # Scale ratio
        belief = gaussian_scaling(src_activity, tgt_activity)
        if is_inhibitory(itype):
            belief = -belief #An activated inhibitory reduces pathway activity.

        node_beliefs[node] = belief

    return float(np.mean(list(node_beliefs.values())))


def process_sample(sample_udp: pd.Series, PATHWAY_GRAPHS):
    """
    Compute pathway activity using graph structures. Sample-specific.
    Gaussian scaling of sum(incoming UDP) and sum(outgoing UDP)
    """
    activities = {}
    
    for pathway, G in PATHWAY_GRAPHS.items():
        graph_activities = []
        for component in nx.weakly_connected_components(G):
            subgraph = G.subgraph(component)
            component_activity = pathway_activity(subgraph, sample_udp)
            graph_activities.append(component_activity)
        
        mean_activity = np.mean(graph_activities)
        activities[pathway] = mean_activity
    
    return activities


def parallel_apply(df, process_sample, PATHWAY_GRAPHS):
    """Applies a function to DataFrame rows in parallel, preserving order."""
    #scaler = MinMaxScaler()
    #df = pd.DataFrame(scaler.fit_transform(df), columns=df.columns, index=df.index)
    n_cores = max(1, mp.cpu_count() - 2) # leave 2 cores free for OS
    process_sample_with_graphs = partial(process_sample, PATHWAY_GRAPHS=PATHWAY_GRAPHS)

    with mp.Pool(n_cores) as pool:
        results = list(
            tqdm(
                pool.imap(process_sample_with_graphs, (row for _, row in df.iterrows())),
                total=len(df),
            )
        )
    return pd.DataFrame(results, index=df.index)


pathway_interactions = parse_pathway_interactions('./data/pathway_relations.csv')


def init_pathway_graphs(PATHWAY_GRAPHS):
    for pathway, interactions in pathway_interactions.items():
        PATHWAY_GRAPHS[pathway] = build_pathway_graph_structure(pathway, interactions)
    print(f"Built {len(PATHWAY_GRAPHS)} pathway graphs")


def calc_activity(udp_file='./data/output_udp.csv', output_file='./data/output_activity.csv'):
    """Main entry: load expression data, run pathway analysis, and save activity matrix."""
    # Initialize graph structures for pathway graphs.
    PATHWAY_GRAPHS = {}
    init_pathway_graphs(PATHWAY_GRAPHS)
    
    udp_df = pd.read_csv(udp_file, sep='\t', index_col=0)
    udp_df.index = udp_df.index.str.lower()
    DEBUG=False
    if DEBUG:
        results_all = []
        for col in udp_df.columns:
            print(f"Processing sample {col}...")
            result = process_sample(udp_df[col], PATHWAY_GRAPHS)
            results_all.append(result)
        results = pd.DataFrame(results_all, index=udp_df.columns).T
    else:
        df_to_process = udp_df.T
        print(f"Processing {len(df_to_process)} samples...")
        results = parallel_apply(df_to_process, process_sample, PATHWAY_GRAPHS).T
    results = results.round(4)
    results.to_csv(output_file)
    print(f"Saved results to {output_file}")
    return results


if __name__ == '__main__':
    calc_activity('./data/TCGACRC_expression-merged.zip')