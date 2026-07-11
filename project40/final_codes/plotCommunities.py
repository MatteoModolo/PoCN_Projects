import os
import pandas as pd
import colorsys
import random
import graph_tool.all as gt


# month to plot
target_month = "2016_01"

network_file = f"final_weighted_networks/as_edges_{target_month}_weighted.tsv.gz"
louvain_file = f"louvain_communities/louvain_nodes_{target_month}.csv"
leiden_file  = f"leiden_communities/leiden_nodes_{target_month}.csv"

output_dir = "bundled_community_plots"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# build graph
# Force pandas to read the ASNs strictly as strings to prevent ".0" float conversions
df = pd.read_csv(network_file, sep='\t', compression='gzip', dtype={'asn1': str, 'asn2': str})

# Strip any accidental whitespace just to be safe
edges = df[['asn1', 'asn2']].apply(lambda x: x.str.strip()).values

g = gt.Graph(directed=False)
asn_prop = g.add_edge_list(edges, hashed=True)

gt.remove_parallel_edges(g)
gt.remove_self_loops(g)

# load other community assignments

try:
    # force pandas to read the CSV ASNs strictly as strings
    louvain_df = pd.read_csv(louvain_file, dtype={'asn': str})
    leiden_df  = pd.read_csv(leiden_file, dtype={'asn': str})
except FileNotFoundError:
    print("Error: Could not find the Louvain or Leiden CSV files for this month.")
    exit()

# Build the dictionaries, ensuring everything is a string
louvain_dict = {str(row['asn']).strip(): int(row['community_id']) for _, row in louvain_df.iterrows()}
leiden_dict  = {str(row['asn']).strip(): int(row['community_id']) for _, row in leiden_df.iterrows()}

louvain_prop = g.new_vertex_property("int")
leiden_prop  = g.new_vertex_property("int")

matched_nodes = 0
for v in g.vertices():
    asn_str = str(asn_prop[v]).strip()
    
    # If a match is found, assign it. Otherwise, put it in community 0.
    l_id = louvain_dict.get(asn_str, 0)
    louvain_prop[v] = l_id
    leiden_prop[v]  = leiden_dict.get(asn_str, 0)
    
    if l_id != 0:
        matched_nodes += 1


def generate_distinct_colors(g, comm_prop):
    color_prop = g.new_vertex_property("vector<double>")
    unique_comms = list(set(comm_prop.get_array()))
    num_comms = len(unique_comms)
    
    random.seed(42)
    random.shuffle(unique_comms)
    
    comm_to_color = {}
    for i, cid in enumerate(unique_comms):
        hue = i / float(num_comms)
        lightness = 0.5 if i % 2 == 0 else 0.7
        saturation = 0.8 if i % 3 == 0 else 1.0
        r, g_c, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        comm_to_color[cid] = (r, g_c, b, 1.0) 
        
    for v in g.vertices():
        color_prop[v] = comm_to_color[comm_prop[v]]
        
    return color_prop

louvain_colors = generate_distinct_colors(g, louvain_prop)
leiden_colors  = generate_distinct_colors(g, leiden_prop)

# louvain plot
louvain_state = gt.BlockState(g, b=louvain_prop)

louvain_output = os.path.join(output_dir, f"louvain_bundled_plot_{target_month}.png")
print(f" -> Rendering Louvain plot to {louvain_output} ")

louvain_state.draw(
    output=louvain_output,
    output_size=(3000, 3000),
    bg_color=[0, 0, 0, 1],
    vertex_size=2.0,
    vertex_fill_color=louvain_colors,
    vertex_color=louvain_colors,
    edge_color=[1, 1, 1, 0.05],
    edge_pen_width=0.3
)

# Leiden plot
leiden_state = gt.BlockState(g, b=leiden_prop)

leiden_output = os.path.join(output_dir, f"leiden_bundled_plot_{target_month}.png")
print(f" -> Rendering Leiden plot to {leiden_output} ")

leiden_state.draw(
    output=leiden_output,
    output_size=(3000, 3000),
    bg_color=[0, 0, 0, 1],
    vertex_size=2.0,
    vertex_fill_color=leiden_colors,
    vertex_color=leiden_colors,
    edge_color=[1, 1, 1, 0.05],
    edge_pen_width=0.3
)

