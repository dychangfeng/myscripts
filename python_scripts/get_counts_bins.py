import numpy as np 
import pandas as pd 

aa = np.linspace(-2000, 2000, 201)
test = [1,1,1,2]
counts, bin_edge=np.histogram(test, bins=aa)
df = pd.DataFrame(np.array(bin_edge))
print(df)
df.to_csv('/Users/Yun/Documents/human_cruciform/bins_edge.txt', sep='\t', index=None)