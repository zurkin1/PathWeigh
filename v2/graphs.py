import os
import sys
import networkx as nx
import plotly.graph_objs as go
import pandas as pd


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(CURRENT_DIR)
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)


# Load and parse interactions into simple pathway_interactions dictionary data structure. e.g.
# pathway_interactions['Adherens junction'][0] = (['baiap2', 'wasf2', 'wasf3', 'wasf1'], 'activation', ['actb', 'actg1'], 'Adherens junction')
def parse_pathway_interactions(relations_file):
    """Parse interactions and assign unique IDs to each"""
    pathway_relations = pd.read_csv(relations_file)
    pathway_relations['source'] = pathway_relations['source'].fillna('').astype(str).str.lower().str.split('*')
    pathway_relations['target'] = pathway_relations['target'].fillna('').astype(str).str.lower().str.split('*')   
    
    interactions_by_pathway = {}
    for idx, row in pathway_relations.iterrows():
        pathway = row['pathway']
        if pathway not in interactions_by_pathway:
            interactions_by_pathway[pathway] = []
        
        # Store interaction with its global ID for fast lookup
        interactions_by_pathway[pathway].append({
            'id': idx,  # Unique interaction ID
            'source': row['source'],
            'type': row['interactiontype'],
            'target': row['target'],
            'pathway': pathway
        })
    
    return interactions_by_pathway


def parse_pathway_interactions(relations_file):
    """
    Parse interactions and assign unique IDs to each
    Load and parse interactions into simple pathway_interactions dictionary data structure. e.g.
    pathway_interactions['Adherens junction'][0] = (['baiap2', 'wasf2', 'wasf3', 'wasf1'], 'activation', ['actb', 'actg1'], 'Adherens junction')
    """
    pathway_relations = pd.read_csv(relations_file)
    pathway_relations['source'] = pathway_relations['source'].fillna('').astype(str).str.lower().str.split('*')
    pathway_relations['target'] = pathway_relations['target'].fillna('').astype(str).str.lower().str.split('*')   
    
    interactions_by_pathway = {}
    for idx, row in pathway_relations.iterrows():
        pathway = row['pathway']
        if pathway not in interactions_by_pathway:
            interactions_by_pathway[pathway] = []
        
        # Store interaction with its global ID for fast lookup
        interactions_by_pathway[pathway].append({
            'id': idx,
            'source': row['source'],
            'type': row['interactiontype'],
            'target': row['target'],
            'pathway': pathway
        })
    
    return interactions_by_pathway


def build_pathway_graph_structure(pathway, interactions):
    """
    Build the static graph structure for a pathway. Stores only topology and gene names no sample-specific data.
    Returns: NetworkX graph with:
    - Nodes: interaction IDs
    - Node attrs: source_genes, target_genes, interaction_type
    - Edge attrs: genes (the shared genes creating this edge)
    - Graph attrs: corridors (list of (source, sink) tuples for top 20 longest paths)
    """
    G = nx.DiGraph()
    
    # Add all nodes first
    for interaction in interactions:
        i_id = interaction['id']
        G.add_node(
            i_id,
            source_genes=interaction['source'],
            target_genes=interaction['target'],
            interaction_type=interaction['type']
        )
    
    # Create edges based on gene sharing (target of i1 → source of i2)
    for int1 in interactions:
        for int2 in interactions:
            if int1['id'] == int2['id']:
                continue
            
            shared_genes = set(int1['target']) & set(int2['source'])
            if shared_genes and shared_genes != {''}:
                # Store all genes that create this connection
                if G.has_edge(int1['id'], int2['id']):
                    # Add to existing gene list
                    G[int1['id']][int2['id']]['genes'].update(shared_genes)
                else:
                    # Create new edge
                    G.add_edge(int1['id'], int2['id'], genes=shared_genes)

    #Calculate SUPERSOURCE and SUPERSINK. Find corridors first.    
    sources = [node for node in G.nodes if G.in_degree(node) == 0]
    sinks   = [node for node in G.nodes if G.out_degree(node) == 0]

    corridors = []
    if sources and sinks:
        # Make temporary acyclic copy for longest path finding.
        G_temp = G.copy()
        while True:
            try:
                cycle = nx.find_cycle(G_temp, orientation='original')
                G_temp.remove_edge(*cycle[0][:2])
            except nx.exception.NetworkXNoCycle:
                break
        
        # Find top 20 longest paths on the temporary DAG
        for _ in range(20):
            try:
                path = nx.dag_longest_path(G_temp)
            except (nx.NetworkXError, nx.NetworkXNotImplemented):
                break
            if len(path) < 2:
                break
            corridors.append((path[0], path[-1]))
            G_temp.remove_edges_from(list(zip(path, path[1:])))
            if G_temp.number_of_edges() == 0:
                break
    G.graph['corridors'] = corridors
    '''
    # Extract unique sources and sinks from corridors
    sources = list(set(src for src, _ in corridors))
    sinks = list(set(snk for _, snk in corridors))
    
    # Add supersource and supersink nodes
    supersource = 'SUPERSOURCE'
    supersink = 'SUPERSINK'
    
    G.add_node(supersource, interaction_type='activation')
    G.add_node(supersink, interaction_type='activation')
    
    # Connect supersource to all corridor sources with infinite capacity
    for source in sources:
        G.add_edge(supersource, source)
    
    # Connect all corridor sinks to supersink with infinite capacity
    for sink in sinks:
        G.add_edge(sink, supersink)
    '''
    #draw_graph(G, pathway)
    #print(f"PATHWAY: {pathway} Corridors: {len(corridors)} edges: {len(G.edges)} nodes: {len(G.nodes)}")
    return G


def draw_graph(G_flow, pathway):
    # Layout for node positions
    pos = nx.spring_layout(G_flow)

    # Nodes
    node_trace = go.Scatter(
        x=[pos[n][0] for n in G_flow.nodes()],
        y=[pos[n][1] for n in G_flow.nodes()],
        text=[str(n) for n in G_flow.nodes()],
        mode='markers+text',
        textposition='top center',
        marker=dict(size=28, color='skyblue', line=dict(width=2, color='black')),
        hoverinfo='text'
    )

    # Edges as lines
    edge_trace = go.Scatter(
        x=[],
        y=[],
        line=dict(width=2, color='gray'),
        mode='lines'
    )

    annotations = []
    for u, v in G_flow.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_trace['x'] += (x0, x1, None)
        edge_trace['y'] += (y0, y1, None)
        # Calculate label position
        label_x = (x0 + x1) / 2
        label_y = (y0 + y1) / 2
        cap = G_flow[u][v]['capacity'] if 'capacity' in G_flow[u][v] else ''
        # Annotate edge label above midpoint
        annotations.append(
            dict(
                x=label_x,
                y=label_y,
                text=str(cap),
                showarrow=False,
                font=dict(size=18),
                yshift=10,
                xanchor='center'
            )
        )
        # Arrow annotation for the edge direction
        annotations.append(
            dict(
                ax=x0, ay=y0,
                x=x1, y=y1,
                xref='x', yref='y',
                axref='x', ayref='y',
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=2,
                arrowcolor='red'
            )
        )

    fig = go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            annotations=annotations
        )
    )
    fig.write_html(f"./data/graphs/{pathway}.html")