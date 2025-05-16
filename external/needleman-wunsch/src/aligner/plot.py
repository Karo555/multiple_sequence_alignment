import matplotlib.pyplot as plt
from typing import List


def plot_matrix(matrix: List[List[int]], path: str, show: bool = False) -> plt.Figure:
    """
    Plot the scoring matrix as a heatmap and save to the specified path.

    :param matrix: 2D list of scores
    :param path: Filepath to save the PNG
    :param show: Whether to display the plot interactively
    :return: Matplotlib Figure object
    """
    fig, ax = plt.subplots()
    cax = ax.imshow(matrix, interpolation="nearest", aspect="auto")
    fig.colorbar(cax, ax=ax)
    ax.set_xlabel("Sequence 2 position")
    ax.set_ylabel("Sequence 1 position")
    ax.set_title("Scoring Matrix Heatmap")
    fig.tight_layout()
    fig.savefig(path)
    if show:
        plt.show()
    plt.close(fig)
    return fig
