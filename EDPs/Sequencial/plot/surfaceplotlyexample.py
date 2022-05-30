import plotly.graph_objects as go
import pandas as pd
import numpy as np

# Read data from a csv
dados = pd.read_csv("bacteria2d_t0.csv")

fig = go.Figure(data=[go.Surface(z=dados,x=dados.x,y=dados.y)])
fig.update_layout(title='Bacteria (t = 0)', autosize=False,
                  width=800, height=600,
                  margin=dict(l=65, r=50, b=65, t=90))
fig.show()

x = np.outer(np.linspace(-2, 2, 30), np.ones(30))
y = x.copy().T
z = np.cos(x ** 2 + y ** 2)

print(x)
print(y)

fig = go.Figure(data=[go.Surface(x=x, y=y, z=z)])
  
fig.show()