import json
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Patch, ConnectionPatch


def load_alphabridge_data(json_path):
    """
    Load AlphaBridge data from the JSON file.

    Args:
        json_path (str): Path to the JSON file.

    Returns:
        dict: Parsed JSON data.
    """
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"File not found: {json_path}")
    
    with open(json_path, "r") as file:
        data = json.load(file)
    return data


def plot_linear_ribbon_diagram_with_interfaces(data, output_path):
    """
    Generate a linear ribbon diagram with interacting interfaces indicated.

    Args:
        data (dict): Parsed AlphaBridge JSON data.
        output_path (str): Path to save the generated plot.
    """
    structures = data.get("structure", [])
    interactions = data.get("interactions_dict", {}).get("interfaces", [])
    if not structures:
        raise ValueError("No structure data found in the JSON file.")

    # Set up the plot
    fig, ax = plt.subplots(figsize=(15, 8))
    ax.set_title("Linear Ribbon Diagram with Interfaces", fontsize=16)
    ax.set_xlabel("Residue Index", fontsize=12)
    ax.set_ylabel("Chains", fontsize=12)

    # Define color ranges based on pLDDT values
    plddt_colors = {
        "Very Low": "#ff7d45",  # Red
        "Low": "#ffdb13",       # Yellow
        "High": "#65cbf3",      # Light Blue
        "Very High": "#0053d6"  # Dark Blue
    }

    # Create a legend for the pLDDT color scheme
    legend_handles = [
        Patch(color=color, label=label) for label, color in plddt_colors.items()
    ]
    ax.legend(handles=legend_handles, loc="upper left", title="Model Confidence")

    # Plot each chain from the structure
    chain_y = 0  # y position for each chain
    chain_positions = {}  # Store positions for chains to plot interactions
    for structure in structures:
        chains = structure.get("chains", {}).get("polymer", [])
        for chain in chains:
            chain_y -= 1  # Move down for the next chain
            residues = chain.get("residues", [])
            chain_label = f"Chain {chain['auth_asym_id']}"

            # Plot residues for this chain
            for residue in residues:
                seq_id = residue["seq_id"]
                plddt = residue["plddt"]

                # Determine color based on pLDDT value
                if plddt < 50:
                    color = plddt_colors["Very Low"]
                elif 50 <= plddt < 70:
                    color = plddt_colors["Low"]
                elif 70 <= plddt < 90:
                    color = plddt_colors["High"]
                else:
                    color = plddt_colors["Very High"]

                # Plot the residue as a rectangle
                ax.plot(
                    [seq_id - 0.5, seq_id + 0.5],  # Residue width
                    [chain_y, chain_y],  # Residue line
                    color=color,
                    lw=5,
                )

            # Add the chain label
            ax.text(-5, chain_y, chain_label, va="center", ha="right", fontsize=10)
            chain_positions[chain["auth_asym_id"]] = chain_y  # Store chain position

    # Highlight interacting interfaces and add connection lines
    for interface in interactions:
        interface_id = interface["interface_id"]
        links = interface.get("links", [])
        for link in links:
            first = link["first"]
            second = link["second"]

            # Get chain and residue range for the first and second interface
            chain_1 = first["asym_id"]
            chain_2 = second["asym_id"]
            range_1 = first["link_range"]
            range_2 = second["link_range"]

            # Highlight interacting ranges in purple
            ax.plot(
                [range_1["start"], range_1["end"]],
                [chain_positions[chain_1], chain_positions[chain_1]],
                color="purple",
                lw=4,
                label=f"Interface {interface_id}" if interface_id not in ax.get_legend_handles_labels()[1] else None,
            )
            ax.plot(
                [range_2["start"], range_2["end"]],
                [chain_positions[chain_2], chain_positions[chain_2]],
                color="purple",
                lw=4,
            )

            # Draw a connection line between the two chains
            con = ConnectionPatch(
                xyA=((range_1["start"] + range_1["end"]) / 2, chain_positions[chain_1]),
                xyB=((range_2["start"] + range_2["end"]) / 2, chain_positions[chain_2]),
                coordsA="data",
                coordsB="data",
                axesA=ax,
                axesB=ax,
                color="gray",
                lw=1,
                alpha=0.7,
            )
            ax.add_artist(con)

    # Adjust axis tick increments for better readability
    ax.xaxis.set_major_locator(MaxNLocator(nbins=20, integer=True))  # Smaller x-axis increments
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))  # Ensure integer ticks for y-axis

    # Finalize the plot
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=chain_y - 1, top=1)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_yticks([])  # Hide y-axis ticks (chains are labeled manually)

    # Save the figure
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Linear Ribbon Diagram with Interfaces saved to {output_path}")


if __name__ == "__main__":
    # Define paths
    json_file = "alphabridge_data.json"  # Replace with your JSON file path
    output_image = "linear_ribbon_diagram_with_interfaces.png"  # Output image path

    try:
        # Load data and generate the plot
        alphabridge_data = load_alphabridge_data(json_file)
        plot_linear_ribbon_diagram_with_interfaces(alphabridge_data, output_image)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")