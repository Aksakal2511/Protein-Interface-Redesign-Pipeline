import os
import re
import plotly.graph_objects as go

# Results directory
file_dir = "${WORKING_DIR}/results"

# Select files
files = [f for f in os.listdir(file_dir) if f.startswith('profile') and f.endswith('.xvg')]

# Read .xvg files
def read_xvg_file(filepath):
    x_data = []
    y_data = []
    with open(filepath, 'r') as file:
        for line in file:
            if not line.startswith(('#', '@')):
                columns = line.split()
                x_data.append(float(columns[0]))
                y_data.append(float(columns[1]))
    return x_data, y_data

# PDB ID data
data_by_pdb = {}

# Function to generate a consistent colour map for residues
def generate_color_map(residues):
    color_palette = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
    ]
    residue_color_map = {}
    for i, residue in enumerate(sorted(residues)):
        residue_color_map[residue] = color_palette[i % len(color_palette)]
    return residue_color_map

# Extract residues from file names
residues = set()
pattern = r"profile_(\w+)_(\w+)_([\d]+)\.xvg"

# Process each file to gather unique residues
for file in files:
    match = re.match(pattern, file)
    if match:
        _, mutation, residue_number = match.groups()
        residues.add(f"{mutation}_{residue_number}")

# Generate a colour map for the residues
residue_color_map = generate_color_map(residues)

# Process each file to group data by PDB ID
for file in files:
    match = re.match(pattern, file)
    if match:
        pdb_id, mutation, residue_number = match.groups()
        label = "{}_{}".format(mutation, residue_number)  # Create a label like 'GLU_456'
        x_data, y_data = read_xvg_file(os.path.join(file_dir, file))
        if pdb_id not in data_by_pdb:
            data_by_pdb[pdb_id] = []
        data_by_pdb[pdb_id].append((label, x_data, y_data))

# Create plots for each PDB
for pdb_id, datasets in data_by_pdb.items():
    fig = go.Figure()
    
    for label, x_data, y_data in datasets:
        fig.add_trace(go.Scatter(
            x=x_data, 
            y=y_data, 
            mode='lines', 
            name=label,
            line=dict(color=residue_color_map[label])  # Use the consistent color for the residue
        ))
    
    fig.update_layout(
        xaxis_title="ξ (nm)",
        xaxis_title_font=dict(size=28, color='black'),
        yaxis_title="E (kcal mol<sup>-1</sup>)",
        yaxis_title_font=dict(size=28, color='black'),  # Double the font size for y-axis label
        xaxis=dict(
            showline=True, 
            linewidth=2, 
            linecolor='black', 
            ticks="outside",
            tickcolor='black',
            tickwidth=2,
            tickfont=dict(size=24, color='black')
        ),
        yaxis=dict(
            showline=True, 
            linewidth=2, 
            linecolor='black', 
            ticks="outside",
            tickcolor='black',
            tickwidth=2,
            tickfont=dict(size=24, color='black')
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        legend=dict(
            font=dict(size=22, color='black'),
            orientation="h",
            x=0.5,
            y=-0.2,
            xanchor="center",
            yanchor="top"
        )
    )
    
    # HTML file
    output_path = os.path.join(file_dir, "{}_profile_plot.html".format(pdb_id))
    fig.write_html(output_path)

    print("Plot saved to {}".format(output_path))

