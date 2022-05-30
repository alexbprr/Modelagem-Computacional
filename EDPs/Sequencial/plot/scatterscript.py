import plotly.express as px
import pandas as pd 
from plotly.subplots import make_subplots
import plotly.graph_objects as go

df = pd.read_csv("bacteria2d_t0.csv")
fig = px.scatter_3d(df, x="x", y="y", z="value", color="value")
fig.show()

df = pd.read_csv("bacteria2d_t1.csv")
fig = px.scatter_3d(df, x="x", y="y", z="value", color="value")
fig.show()

df = pd.read_csv("bacteria2d_t5.csv")
fig = px.scatter_3d(df, x="x", y="y", z="value", color="value")
fig.show()

df = pd.read_csv("bacteria2d_t10.csv")
fig = px.scatter_3d(df, x="x", y="y", z="value", color="value")
fig.update_layout(title="Bacteria (t=10)")
fig.show()

df = pd.read_csv("bacteria2d_t50.csv")
fig = px.scatter_3d(df, x="x", y="y", z="value", color="value")
fig.update_layout(title="Bacteria (t=50)")
fig.show()

