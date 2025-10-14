import plotly.graph_objects as go

def plot_gc_content(positions, gc_values):
    """Creates a line plot of GC content vs. sequence position."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=positions, y=gc_values, mode='lines', name='GC Content'))
    
    fig.update_layout(
        title="GC Content Sliding Window Analysis",
        xaxis_title="Position (bp)",
        yaxis_title="GC Content (%)",
        yaxis=dict(range=[0, 100])
    )
    return fig
import pandas as pd

def plot_codon_heatmap(codon_dict):
    """Creates a heatmap of codon usage."""
    if not codon_dict:
        return go.Figure()

    # Convert dictionary to a pandas DataFrame for easier plotting
    df = pd.DataFrame(list(codon_dict.items()), columns=['Codon', 'Count'])
    
    # Create a pivot table for the heatmap
    codons = sorted(df['Codon'].unique())
    bases = ['A', 'C', 'G', 'T']
    heatmap_df = pd.DataFrame(index=[b1+b2 for b1 in bases for b2 in bases], columns=bases)

    for codon, count in codon_dict.items():
        if len(codon) == 3:
            row, col = codon[:2], codon[2]
            heatmap_df.loc[row, col] = count

    heatmap_df = heatmap_df.fillna(0) # Fill NaNs with 0

    fig = go.Figure(data=go.Heatmap(
                   z=heatmap_df.values,
                   x=heatmap_df.columns,
                   y=heatmap_df.index,
                   hoverongaps=False,
                   colorscale='Viridis'))
    
    fig.update_layout(title="Codon Usage Frequency Heatmap")
    return fig
def plot_kmer_distribution(kmer_dict):
    """Creates a bar chart of the top 10 most frequent k-mers."""
    if not kmer_dict:
        return go.Figure()

    # Convert dictionary to a DataFrame and get the top 10
    df = pd.DataFrame(list(kmer_dict.items()), columns=['k-mer', 'Count'])
    df = df.sort_values(by='Count', ascending=False).head(10)

    fig = go.Figure([go.Bar(x=df['k-mer'], y=df['Count'])])
    fig.update_layout(
        title="Top 10 Most Frequent k-mers",
        xaxis_title="k-mer",
        yaxis_title="Frequency Count"
    )
    return fig