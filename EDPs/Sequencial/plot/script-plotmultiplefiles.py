import glob
import numpy as np
import pandas as pd
from typing import Iterable
import matplotlib.pyplot as plt
import sys

def read_np_from_csv(filename: str) -> 'tuple[np.array, ...]':
    return pd.read_csv(filename).to_numpy().T

def write_np_to_csv(filename: str, *data, columns: 'Iterable[str]'=None) -> None:
    data_matrix = np.matrix(data).T
    df = pd.DataFrame(data_matrix, columns=columns)
    df.to_csv(filename, index=False)

filename = sys.argv[1]
results = glob.glob(filename + "_t*.csv")
print(results)
for ind, result in enumerate(results):

    # Supondo que o CSV está da seguinte forma:
    # valor_x, valor_y, *multiplos_outros_valores, ...
    #x, values = read_np_from_csv(result)

    df = pd.read_csv(result)

    fig = plt.figure()
    # Plot x line and y
    plt.plot(df["x"], df["value"], label=ind)
    fig.savefig(result + ".svg", format="svg", bbox_inches="tight")
    #plt.show()