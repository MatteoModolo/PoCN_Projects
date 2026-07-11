import os
import glob
import time
import re
import pandas as pd
import graph_tool.all as gt


# data and dir
input_dir = "final_weighted_networks"
output_summary_csv = "dcsbm_mdl_summary.csv" # Updated name
comm_output_dir = "dcsbm_mdl_blocks"         # Updated name

target_plot_month = "2016_01"

if not os.path.exists(comm_output_dir):
    os.makedirs(comm_output_dir)

network_files = sorted(glob.glob(os.path.join(input_dir, "as_edges_*.tsv.gz")))

if len(network_files) == 0:
    print(f"No network files found in '{input_dir}'.")
    exit()

results_list = []


for file_path in network_files:
    
    filename = os.path.basename(file_path)
    match = re.search(r"(\d{4}_\d{2})", filename)
    if not match:
        continue
    yyyy_mm = match.group(1)
    
    print(f"\n[{yyyy_mm}] Loading network...")
    df = pd.read_csv(file_path, sep='\t', compression='gzip')
    
    # Build the graph-tool Graph
    edges = df[['asn1', 'asn2']].astype(str).values
    
    g = gt.Graph(directed=False)
    asn_prop = g.add_edge_list(edges, hashed=True)
    
    gt.remove_parallel_edges(g)
    gt.remove_self_loops(g)
    
    # Run MDL optimization with degree correction
    start_time = time.time()
    
    # state_args=dict(deg_corr=True) triggers the Degree-Corrected SBM
    state = gt.minimize_blockmodel_dl(g, state_args=dict(deg_corr=True))
    
    calc_time = round((time.time() - start_time) / 60, 2)
    
    # Extract the Metrics
    mdl_score = state.entropy()
    num_blocks = state.get_B()
    num_nodes = g.num_vertices()
    
    print(f"[{yyyy_mm}]| Nodes: {num_nodes} | Blocks: {num_blocks} | MDL (Entropy): {mdl_score:.1f} | Time: {calc_time} mins")

    block_map = state.get_blocks()
    
    asns = [asn_prop[v] for v in g.vertices()]
    blocks = [block_map[v] for v in g.vertices()]
    
    comm_df = pd.DataFrame({"asn": asns, "dcsbm_block_id": blocks})
    comm_df.to_csv(os.path.join(comm_output_dir, f"dcsbm_nodes_{yyyy_mm}.csv"), index=False)
    
    # visulization
    if yyyy_mm == target_plot_month:
        print(f"   -> [Visualization] Generating native DC-SBM edge-bundled plot for {yyyy_mm}...")
        plot_filename = os.path.join(comm_output_dir, f"dcsbm_mdl_plot_{yyyy_mm}.png")
        
        state.draw(
            output=plot_filename,
            output_size=(3000, 3000),
            bg_color=[1, 1, 1, 1],
            vertex_size=2,
            edge_pen_width=0.5
        )
        print(f"   -> [Visualization] DC-SBM Block plot saved to {plot_filename}")

    # Store metrics
    results_list.append({
        "month": yyyy_mm,
        "total_nodes": num_nodes,
        "dcsbm_blocks": num_blocks,
        "mdl_entropy_score": mdl_score,
        "compute_time_mins": calc_time
    })

# save data
final_summary = pd.DataFrame(results_list)
final_summary.to_csv(output_summary_csv, index=False)


print(f" Master summary saved to: {output_summary_csv}")
print(f" Individual block assignments saved in: {comm_output_dir}/")
