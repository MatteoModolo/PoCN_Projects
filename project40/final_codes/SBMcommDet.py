import os
import glob
import time
import re
import pandas as pd
import graph_tool.all as gt

# Directories and files
input_dir = "final_weighted_networks"
output_summary_csv = "sbm_mdl_summary.csv"
comm_output_dir = "sbm_mdl_blocks"

target_plot_month = "2016_01"

if not os.path.exists(comm_output_dir):
    os.makedirs(comm_output_dir)

network_files = sorted(glob.glob(os.path.join(input_dir, "as_edges_*.tsv.gz")))

if len(network_files) == 0:
    print(f"No network files found in '{input_dir}'.")
    exit()

results_list = []


# MDL Stochastic Block Model analysis
for file_path in network_files:
    
    filename = os.path.basename(file_path)
    match = re.search(r"(\d{4}_\d{2})", filename)
    if not match:
        continue
    yyyy_mm = match.group(1)
    
    df = pd.read_csv(file_path, sep='\t', compression='gzip')

    #create graph
    edges = df[['asn1', 'asn2']].astype(str).values
    
    g = gt.Graph(directed=False)
    # hashed=True safely maps ASNs to internal integer vertex IDs
    asn_prop = g.add_edge_list(edges, hashed=True)
    
    
    # Run the Minimum Description Length optmization
    #take time to include even this as a tradeoff
    start_time = time.time()
     
    # It automatically finds the optimal number of blocks minimizing descriprion length
    state = gt.minimize_blockmodel_dl(g)
    
    calc_time = round((time.time() - start_time) / 60, 2)
    
    # Extract the metrics
    mdl_score = state.entropy()
    num_blocks = state.get_B()
    num_nodes = g.num_vertices()
    
    print(f"[{yyyy_mm}] | Nodes: {num_nodes} | Blocks: {num_blocks} | MDL (Entropy): {mdl_score:.1f} | Time: {calc_time} mins")
    
    # Extract assignments
    # state.get_blocks() returns a property map of the block assigned to each vertex
    block_map = state.get_blocks()
    
    asns = [asn_prop[v] for v in g.vertices()]
    blocks = [block_map[v] for v in g.vertices()]
    
    comm_df = pd.DataFrame({"asn": asns, "sbm_block_id": blocks})
    comm_df.to_csv(os.path.join(comm_output_dir, f"sbm_nodes_{yyyy_mm}.csv"), index=False)
    
    # plots
    if yyyy_mm == target_plot_month:
        print(f"Generating native SBM edge-bundled plot for {yyyy_mm}...")
        plot_filename = os.path.join(comm_output_dir, f"sbm_mdl_plot_{yyyy_mm}.png")
        
        # use graphtool native visualization (due to this the layouts of the other methods are adapted to this)
        state.draw(
            output=plot_filename,
            output_size=(3000, 3000),
            bg_color=[1, 1, 1, 1],
            vertex_size=2,
            edge_pen_width=0.5
        )
        print(f"MDL Block plot saved to {plot_filename}")

    # 5. Store metrics
    results_list.append({
        "month": yyyy_mm,
        "total_nodes": num_nodes,
        "sbm_blocks": num_blocks,
        "mdl_entropy_score": mdl_score,
        "compute_time_mins": calc_time
    })

# save data
final_summary = pd.DataFrame(results_list)
final_summary.to_csv(output_summary_csv, index=False)

print(f" Master summary saved to: {output_summary_csv}")
print(f" Individual block assignments saved in: {comm_output_dir}/")
