import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc
from scipy.sparse import find
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge
from scipy.sparse import find
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()

# Function to create a pie chart at a specific location
def _pie_chart(ax, ratios, center, size, colors):
    ratios = ratios / ratios.sum()  # Ensure the ratios sum to 1
    angle = 0
    for r, color in zip(ratios, colors):
        wedge = Wedge(center, size, angle, angle + 360 * r, facecolor=color)
        ax.add_patch(wedge)
        angle += 360 * r

class link_comm:
    """
    Class to run link communities on a graph of single-cell AnnData object and visualize the results.
    """
    def __init__(self, adata:sc.AnnData, graph:str="connectivities"):
        """
        Initialize the link_comm object with an AnnData object and a graph.

        Args:
            adata (AnnData): Input single-cell AnnData object
            graph (str): Name of the graph in adata.obsp to use for link communities detection. Defaults to "connectivities".
        """
        adjacency_matrix = adata.obsp[graph]

        # Convert adjacency matrix to edge list
        rows, cols, values = find(adjacency_matrix)
        edge_list = list(zip(rows, cols, values))


        # Convert edge list to dataframe
        edge_df = pd.DataFrame(edge_list, columns=["node1", "node2", "weight"])
        obs_mapping = dict(zip(range(adata.n_obs), adata.obs_names))
        for col in ["node1", "node2"]:
            edge_df[col] = edge_df[col].map(obs_mapping)
        self.edge_df = edge_df
        self.proportions_df = None
    
    def run_link_comm(self, find_meta_communities=True, hcmethod:str="average", use_all_edges=False, check_duplicates=True, edglim=10000):
        """
        Run link communities on the graph.\n
        
        Once the link communities are detected, the results are stored in the object with the following attributes:
        - link_comm_result: A dictionary with the following keys:
            - node_clusters: A pandas DataFrame with the node clusters
            - link_clusters: A pandas DataFrame with the link clusters\n
            
        If find_meta_communities is set to True, the meta-communities are also stored in the object with the following attributes:
        - meta_communities: A dictionary with the following keys:
            - node_clusters: A pandas DataFrame with the node clusters
            - link_clusters: A pandas DataFrame with the link clusters

        Args:
            find_meta_communities (bool, optional): Whether to find meta-communities. Defaults to True.
            hcmethod (str, optional): Hierarchical clustering method. Defaults to "average".
            use_all_edges (bool, optional): Whether to use all edges. Defaults to False.
            check_duplicates (bool, optional): Whether to check for duplicates. Defaults to True.
            edglim (int, optional): An integer value indicating the largest number of edges permissible for the hierarchical clustering to be handled in memory. Defaults to 10000.
        """

        # import linkcomm R packages
        linkcomm = importr("linkcomm")
        
        # Run linkcomm
        lc = linkcomm.getLinkCommunities(pandas2ri.py2rpy(self.edge_df), hcmethod, use_all_edges, edglim, check_duplicates=check_duplicates, plot=False)
        
        # Extract link communities
        self.link_comm_result = {}
        self.link_comm_result['node_clusters'] = ro.conversion.rpy2py(lc[4])
        self.link_comm_result['link_clusters'] = self.edge_df.copy()
        link_cluster_map = {i: int(node+1) for node, cluster in enumerate(lc[5]) for i in cluster.astype(int)}
        self.link_comm_result['link_clusters']['cluster'] = self.edge_df.index.map(link_cluster_map)
        
        # Extract meta-communities
        if find_meta_communities:
            mc = linkcomm.meta_communities(lc)
            self.meta_communities = {}
            self.meta_communities['node_clusters'] = ro.conversion.rpy2py(mc[4])
            self.meta_communities['link_clusters'] = self.edge_df.copy()
            link_cluster_map = {i: int(node+1) for node, cluster in enumerate(mc[5]) for i in cluster.astype(int)}
            self.meta_communities['link_clusters']['cluster'] = self.edge_df.index.map(link_cluster_map)
            
    def plot_link_comm(self, adata:sc.AnnData, use_meta:bool=True, embedding:str='X_umap', figsize=(5, 5)):
        """Plot the link communities on a 2D-embedding plot of the single-cell data.

        Args:
            adata (AnnData): Input single-cell AnnData object
            use_meta (bool, optional): Whether to use meta-communities for plotting. Defaults to True.
            embedding (str, optional): Key for the embedding stored in adata.obsm. Defaults to 'X_umap'.
            figsize (tuple, optional): Figure size. Defaults to (5, 5).

        Returns:
            tuple: A tuple containing the figure and axis objects.
        """
        
        # select the edge dataframe
        edge_df = self.meta_communities['link_clusters'].copy() if use_meta else self.link_comm_result['link_clusters'].copy()
        edge_df = edge_df.dropna()
        
        # If proportions_df is not available, calculate it
        if self.proportions_df is None:
            # Aggregate link cluster information for each node
            node_clusters = defaultdict(list)

            # extract node clusters
            for _, row in edge_df.iterrows():
                node1, node2, cluster = row['node1'], row['node2'], row['cluster']
                node_clusters[node1].append(cluster)
                node_clusters[node2].append(cluster)
                
            # Calculate proportion of each cluster for each node
            node_cluster_proportions = {}
            for node, clusters in node_clusters.items():
                total = len(clusters)
                proportions = pd.Series(clusters).value_counts() / total
                node_cluster_proportions[node] = proportions.to_dict()
            # Convert to DataFrame
            proportions_df = pd.DataFrame(node_cluster_proportions).T.fillna(0)
            proportions_df.columns = proportions_df.columns.astype(int)
            proportions_df = proportions_df[proportions_df.columns.sort_values()]
            self.proportions_df = proportions_df.copy()
        else:
            proportions_df = self.proportions_df.copy()
        
        # UMAP coordinates
        X_coords = adata.obsm[embedding]

        # Colors for the clusters
        colors = plt.cm.tab20b(np.linspace(0, 1, proportions_df.shape[1]))

        fig, ax = plt.subplots(figsize=figsize)
        fig.gca().set_aspect('equal', adjustable='box')

        # Scatter plot for UMAP coordinates
        ax.scatter(X_coords[:, 0], X_coords[:, 1], s=10, facecolors='none', edgecolors='k', alpha=0.5, lw=0.5)

        # Add pie charts to the plot
        for (i, (x, y)) in enumerate(X_coords):
            node = adata.obs.index[i]
            if node in proportions_df.index:
                ratios = proportions_df.loc[node].values
                _pie_chart(ax, ratios, (x, y), size=0.15, colors=colors)
                
        # Plot edges with color based on clusters
        edge_df_node_idx = edge_df.copy()
        for col in ['node1', 'node2']:
            edge_df_node_idx[col] = edge_df_node_idx[col].map(dict(zip(adata.obs_names, range(adata.n_obs))))
        edge_colors = edge_df['cluster']  # Cluster labels
        unique_clusters = np.unique(edge_colors)
        for cluster in unique_clusters:
            cluster_edges = edge_df_node_idx.dropna()[edge_colors == cluster]
            for i, j, in zip(cluster_edges['node1'], cluster_edges['node2']):
                ax.plot([X_coords[int(i), 0], X_coords[int(j), 0]],
                        [X_coords[int(i), 1], X_coords[int(j), 1]],
                        c=colors[int(cluster-1)], alpha=0.3, lw=0.5, zorder=-1)

        # Add legend
        legend_patches = [mpatches.Patch(color=colors[i], label=f'{c}') for i,c in enumerate(proportions_df.columns)]
        ax.legend(handles=legend_patches, loc='center left',frameon=False,
                  bbox_to_anchor=(1, 0.5), ncol=2, fontsize=8, columnspacing=0.5)

        # Set plot title and axis labels
        ax.set_title('Link_Cluster')
        ax.set_xlabel(f'{embedding[2:].upper()}1')
        ax.set_ylabel(f'{embedding[2:].upper()}2')
        ax.axis('off')
        ax.grid(False)
        
        return fig, ax