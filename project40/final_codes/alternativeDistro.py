import powerlaw
import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt          
import matplotlib.dates as mdates        

input_dir = "degree_distribution"
stats_file = "network_monthly_statistics.csv"
output_csv = "ultimate_distribution_confrontation.csv"

# Load the network statistics to reconstruct the raw data
try:
    stats_df = pd.read_csv(stats_file)
except FileNotFoundError:
    print(f"Error: Could not find {stats_file}")
    exit()

degree_files = glob.glob(os.path.join(input_dir, "degree_dist_*.csv"))
if not degree_files:
    print("No degree distribution files found ")
    exit()

results_list = []

# loop on networks and fit all models
for filepath in degree_files:
    filename = os.path.basename(filepath)
    yyyy_mm = filename.replace("degree_dist_", "").replace(".csv", "")
    
    # read data
    deg_df = pd.read_csv(filepath)
    
    # take total nodes
    month_total_nodes = stats_df.loc[stats_df['month'] == yyyy_mm, 'unique_nodes'].values[0]
        
    # Multiply probability by total nodes and reconstruct the array
    counts = np.round(deg_df['probability'] * month_total_nodes).astype(int)
    raw_degrees = np.repeat(deg_df['degree'].values, counts)
    
    # Filter out 0
    valid_degrees = raw_degrees[raw_degrees > 0]
    
    # Fit of all distro and comparison
    fit = powerlaw.Fit(valid_degrees, discrete=True, verbose=False)
    
    R_ln, p_ln = fit.distribution_compare('power_law', 'lognormal', normalized_ratio=True)
    R_exp, p_exp = fit.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    R_trunc, p_trunc = fit.distribution_compare('power_law', 'truncated_power_law', normalized_ratio=True)
    R_str, p_str = fit.distribution_compare('power_law', 'stretched_exponential', normalized_ratio=True)
    
    # Store results
    results_list.append({
        'month': yyyy_mm,
        'x_min': fit.xmin,
        'alpha': fit.power_law.alpha,
        
        # Log-Normal
        'LN_mu': fit.lognormal.mu,
        'LN_sigma': fit.lognormal.sigma,
        'R_PL_vs_LN': R_ln,
        'p_PL_vs_LN': p_ln,
        
        # Exponential
        'EXP_lambda': fit.exponential.Lambda,
        'R_PL_vs_EXP': R_exp,
        'p_PL_vs_EXP': p_exp,
        
        # Powelaw with exponential cutoff 
        'TRUNC_alpha': fit.truncated_power_law.alpha,
        'TRUNC_lambda': fit.truncated_power_law.Lambda,
        'R_PL_vs_TRUNC': R_trunc,
        'p_PL_vs_TRUNC': p_trunc,
        
        # Stretched exp
        'STR_beta': fit.stretched_exponential.beta,
        'STR_lambda': fit.stretched_exponential.Lambda,
        'R_PL_vs_STR': R_str,
        'p_PL_vs_STR': p_str
    })
    
    print(f" PL vs Cutoff: R={R_trunc:.2f} (p={p_trunc:.3f})")

# save data 
final_df = pd.DataFrame(results_list)

# Mathematically sort the dataframe chronologically to ensure the line graph flows left to right
final_df = final_df.sort_values('month')
final_df.to_csv(output_csv, index=False)
print(f"\nData extraction complete. Saved to '{output_csv}'")


# Convert the YYYY_MM strings into actual datetime objects for a perfectly smooth X-axis
final_df['date'] = pd.to_datetime(final_df['month'], format='%Y_%m')

plt.figure(figsize=(10, 6))

# Plot the 4 Log-Likelihood Ratios with distinct, colorblind-friendly colors
plt.plot(final_df['date'], final_df['R_PL_vs_LN'], label='vs Lognormal', color='#377EB8', lw=2)
plt.plot(final_df['date'], final_df['R_PL_vs_EXP'], label='vs Exponential', color='#E41A1C', lw=2)
plt.plot(final_df['date'], final_df['R_PL_vs_TRUNC'], label='vs Truncated PL', color='#4DAF4A', lw=2)
plt.plot(final_df['date'], final_df['R_PL_vs_STR'], label='vs Stretched Exp', color='#984EA3', lw=2)

# The most important line on the chart: the Zero Line (Tie)
plt.axhline(0, color='black', linestyle='--', linewidth=1.5, zorder=5)

# Formatting and aesthetics
plt.title("Log-Likelihood Ratio (R) of Power-Law vs Alternative Distributions", fontsize=14, fontweight='bold')
plt.ylabel("Log-Likelihood Ratio (R)", fontsize=12)
plt.xlabel("Timeframe", fontsize=12)

# Clean up the X-axis to only show every 2 years (base=2), preventing messy text overlap
plt.gca().xaxis.set_major_locator(mdates.YearLocator(base=2))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

# Add a clean legend
plt.legend(title="Power-Law Matchups", loc='lower left', framealpha=0.9)
plt.grid(True, axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

# Save the plot
plot_filename = "llr_comparison_timeline.png"
plt.savefig(plot_filename, dpi=300, facecolor='white')
print(f"Done! Clean timeline plot saved as '{plot_filename}'")
