import plotly.express as px
import pandas as pd 
from plotly.subplots import make_subplots
import plotly.graph_objects as go

for t in [0,1,5,10,20,30,40,50]:
    df = pd.read_csv("bacteria2d_t" + str(t) + ".csv")
    fig = px.scatter_3d(df, x="x", y="y", z="value", color="value")
    fig.update_layout(title="Bacteria (t=" +str(t) +")")
    fig.show()
    #fig.write_image("images/bacteria2d_t"+ str(t) + ".svg") 

