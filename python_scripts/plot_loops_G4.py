

import pandas as pd
hg38 = pd.read_csv('/uufs/chpc.utah.edu/common/home/u0674686/ncbi_genomes/hg38_all_g4_loops.csv')

hg38.head()

hg38 = hg38[['g_runs', 'loop1', 'loop2', 'loop3', 'loop4', 'loop5']] ## remove the unwanted columns

hg38.describe()

hg38_g_runs = list(hg38.g_runs)

len(hg38.g_runs)

hg38_g_runs_unique = list(set(hg38.g_runs)) ## unique 30391

hg38_g_runs_unique = hg38_g_runs_unique.sort()

hg38_g_runs_unique = sort(hg38_g_runs_unique)




hg38_loops_group = hg38.groupby('g_runs')

help(hg38_loops_group)

help(hg38_loops_group.size())

g_runs_dist = hg38_loops_group.size().sort_values(ascending = False)## len(g_runs_dist)=30391

len(g_runs_dist)

import seaborn as sns
%matplotlib inline




g_runs_dist = g_runs_dist.to_frame() ## convert pandas series to a data frame
g_runs_dist['G_runs'] = g_runs_dist.index
#g_runs_dist = g_runs_dist.reset_index()
g_runs_dist.columns = ['Counts', 'G_runs']


g_runs_dist['g_runs_size'] = [len(x)/2 for x in g_runs_dist.G_runs] # add the how many runs of G as a column 

g_runs_dist

distribution ={}
for i in range(4,20):
    g_n=g_runs_dist[g_runs_dist.g_runs_size>=i].Counts.sum()
    distribution[i]=float(g_n*100/687761)
print distribution
distribution=pd.DataFrame.from_dict(distribution, orient='index')
distribution = distribution.reset_index()
distribution.columns=['G_runs', 'Percentage']

import seaborn as sns
sns.set_style('dark') ## set themes
sns.set_style('ticks')
fig, ax = plt.subplots()
# the size of A4 paper
fig.set_size_inches(8, 8)
plt.xticks(rotation=90)
sns.set_style('ticks')

#sns.set_context('talk', rc={"font.size":15,"axes.titlesize":15,"axes.labelsize":15})
#sns.set(font_scale=1)
#colors = sns.color_palette("husl", n_colors=2)
n_g_runs = sns.barplot(data = distribution, x = 'G_runs', y = 'Percentage')
n_g_runs.set_xlabel("Number of G runs", fontsize=20)
n_g_runs.set_ylabel("Percentage of G4s in all G4s", fontsize = 20)
n_g_runs.tick_params(labelsize=16)
n_g_runs.axes.set_title("The number of G runs in all G4s in human genome", fontsize=24)
plt.savefig("/uufs/chpc.utah.edu/common/home/u0674686/Desktop/number_g_runs_hg38.png", dpi =500)
n_g_runs

g_runs_dist.sort_values(by=['g_runs_size'], ascending=False)

current_palette = sns.color_palette()
sns.palplot(current_palette)


sns.set_style('dark') ## set themes
sns.set_style('ticks')
fig, ax = plt.subplots()
# the size of A4 paper
fig.set_size_inches(10, 8)
plt.xticks(rotation=60)
sns.set_style('ticks')
#sns.set_context('talk', rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":5})
#sns.set(font_scale=1)
colors = sns.color_palette("husl", n_colors=2) # Available seaborn palette names: deep, muted, bright, pastel, dark, colorblind, hls, husl, any named matplotlib palette, list of colors
color_dict = dict(zip(unique(g_runs_dist.g_runs_size[0:20]), colors))
#palette_dir = {1:'blue', 2:'red', 3:'purple', 4:'green'}
palette = (g_runs_dist[0:20])['g_runs_size'].map(color_dict) 
x1=sns.barplot(data = g_runs_dist[0:20], x = 'G_runs', y = 'Counts', palette = palette)
x1.set_xlabel("G4 patterns", fontsize=20)
x1.set_ylabel("The number of G4s", fontsize = 20)
x1.tick_params(labelsize=16)
x1.axes.set_title("The number of G4s in each specific pattern", fontsize=24)
plt.tight_layout()
plt.savefig("/uufs/chpc.utah.edu/common/home/u0674686/Desktop/hg38_g_runs_distribution_40.png", dpi =500)

hg38_33333.describe()

hg38_33333 = (hg38[hg38.g_runs == '3:3:3:3:3:'])[['g_runs', 'loop1', 'loop2', 'loop3', 'loop4']]
hg38_3333 = (hg38[hg38.g_runs == '3:3:3:3:'])[['g_runs', 'loop1', 'loop2', 'loop3']]

hg38_33333.describe()

hg38_3333.describe()

hg38_33333.loop2.value_counts()

melt_hg38_3333 = pd.melt(hg38_3333, id_vars = ['g_runs'], var_name = 'all_loop') ## melt all three loops to one column

melt_hg38_33333 = pd.melt(hg38_33333, id_vars = ['g_runs'], var_name = 'all_loop') ## melt all three loops to one column

all_loop_hg38_33333=(melt_hg38_33333.value.value_counts()).to_frame()
all_loop_hg38_33333['loop']=all_loop_hg38_33333.index

all_loop_hg38_33333['loop_length'] = [len(x) for x in all_loop_hg38_33333.loop] ## add loop length to dataframe

all_loop_hg38_33333.head()

all_loop_hg38_3333 = (melt_hg38_3333.value.value_counts()).to_frame()

all_loop_hg38_3333['loop'] = all_loop_hg38_3333.index

[len(x) for x in all_loop_hg38_3333.loop]

all_loop_hg38_3333['loop_length'] = [len(x) for x in all_loop_hg38_3333.loop]

all_loop_hg38_3333.loop[0:30]



sort(all_loop_hg38_3333.loop_length[0:20])

sns.set_style('dark') ## set themes
sns.set_style('ticks')
fig, ax = plt.subplots()
# the size of A4 paper
fig.set_size_inches(10, 8)
plt.xticks(rotation=60)
sns.set_style('ticks')
sns.set_context('talk')
plt.title('hg38_3333_loop_size')
sns.set(font_scale=1)
colors = sns.color_palette("husl", n_colors=len(set(all_loop_hg38_3333.loop_length[0:30]))) # Available seaborn palette names: deep, muted, bright, pastel, dark, colorblind, hls, husl, any named matplotlib palette, list of colors
color_dict = dict(zip(unique(all_loop_hg38_3333.loop_length[0:30]), colors))
#palette_dir = {1:'blue', 2:'red', 3:'purple', 4:'green'}
palette = (all_loop_hg38_3333[0:30])['loop_length'].map(color_dict) # ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
g = sns.barplot(data = all_loop_hg38_3333[0:30], x = 'loop', y = 'value', palette = palette)
g.set_xlabel("Loops", fontsize=20)
g.set_ylabel("Counts of each specific loop", fontsize = 20)
g.tick_params(labelsize=16)
g.axes.set_title("The distribution of different loops in 3333 G4s ", fontsize=24)
plt.tight_layout()
plt.savefig('/uufs/chpc.utah.edu/common/home/u0674686/Desktop/hg38_3333_loop_distribution', dpi = 500)

unique(all_loop_hg38_33333.loop_length[0:30])

import seaborn as sns
sns.set_style('dark') ## set themes
sns.set_style('ticks')
fig, ax = plt.subplots()
# the size of A4 paper
fig.set_size_inches(10, 8)
plt.xticks(rotation=90)
sns.set_style('ticks')
sns.set_context('talk')
plt.title('hg38_33333_loops')
sns.set(font_scale=1)
colors = sns.color_palette("husl", n_colors=len(zip(unique(all_loop_hg38_33333.loop_length[0:30])))) # Available seaborn palette names: deep, muted, bright, pastel, dark, colorblind, hls, husl, any named matplotlib palette, list of colors
color_dict = dict(zip(unique(all_loop_hg38_33333.loop_length[0:30]), colors))
#palette_dir = {1:'blue', 2:'red', 3:'purple', 4:'green'}
palette = (all_loop_hg38_33333[0:30])['loop_length'].map(color_dict) # ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
g = sns.barplot(data = all_loop_hg38_33333[0:30], x = 'loop', y = 'value', palette = palette)
g.set_xlabel("Loops", fontsize=20)
g.set_ylabel("Counts of each specific loop", fontsize = 20)
g.tick_params(labelsize=16)
g.axes.set_title("The distribution of different loops in 33333 G4s ", fontsize=24)
plt.tight_layout()
plt.savefig('/uufs/chpc.utah.edu/common/home/u0674686/Desktop/hg38_33333_loop_distribution', dpi = 500)

hg38_4333 = (hg38[hg38.g_runs == '4:3:3:3:'])[['g_runs', 'loop1', 'loop2', 'loop3']]
melt_hg38_4333 = pd.melt(hg38_4333, id_vars = ['g_runs'], var_name = 'all_loop') ## melt all three loops to one column
all_loop_hg38_4333=(melt_hg38_4333.value.value_counts()).to_frame()
all_loop_hg38_4333['loop']=all_loop_hg38_4333.index
all_loop_hg38_4333['loop_length'] = [len(x) for x in all_loop_hg38_4333.loop]

help(g.get_legend_handles_labels())

color_dict

sns.set_style('dark') ## set themes
sns.set_style('ticks')
fig, ax = plt.subplots()
# the size of A4 paper
fig.set_size_inches(10, 8)
plt.xticks(rotation=90)
sns.set_style('ticks')
sns.set_context('talk')
plt.title('hg38_4333_loops')
sns.set(font_scale=1)
colors = sns.color_palette("husl", n_colors=len(unique(all_loop_hg38_4333.loop_length[0:30]))) # Available seaborn palette names: deep, muted, bright, pastel, dark, colorblind, hls, husl, any named matplotlib palette, list of colors
color_dict = dict(zip(unique(all_loop_hg38_4333.loop_length[0:30]), colors))
palette = (all_loop_hg38_4333[0:30])['loop_length'].map(color_dict) # ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
g = sns.barplot(data = all_loop_hg38_4333[0:30], x = 'loop', y = 'value', palette = palette)
#g=sns.pointplot(data = all_loop_hg38_4333[0:20], x='loop', y='value', hue='loop_length')
g.set_xlabel("Loops", fontsize=20)
g.set_ylabel("Counts of each specific loop", fontsize = 20)
g.tick_params(labelsize=16)
g.axes.set_title("The distribution of different loops in 4333 G4s ", fontsize=24)
handles, labels = g.get_legend_handles_labels()
g.legend()
plt.tight_layout()
plt.savefig('/uufs/chpc.utah.edu/common/home/u0674686/Desktop/hg38_4333_loop_distribution', dpi = 500)

hg38.head()

hg38_3loops = hg38[['g_runs', 'loop1', 'loop2', 'loop3']]
melt_hg38_3loops = pd.melt(hg38_3loops, id_vars = ['g_runs'], var_name = 'all_loop') ## melt all three loops to one column
all_loop_hg38_3loops=(melt_hg38_3loops.value.value_counts()).to_frame()
all_loop_hg38_3loops['loop']=all_loop_hg38_3loops.index
all_loop_hg38_3loops['loop_length'] = [len(x) for x in all_loop_hg38_3loops.loop]

all_loop_hg38_3loops.head()

sns.set_style('dark') ## set themes
sns.set_style('ticks')
fig, ax = plt.subplots()
# the size of A4 paper
fig.set_size_inches(10, 8)
plt.xticks(rotation=60)
sns.set_style('ticks')
sns.set_context('talk')
colors = sns.color_palette("husl", n_colors=len(unique(all_loop_hg38_3loops.loop_length[0:20]))) # Available seaborn palette names: deep, muted, bright, pastel, dark, colorblind, hls, husl, any named matplotlib palette, list of colors
color_dict = dict(zip(unique(all_loop_hg38_3loops.loop_length[0:20]), colors))
palette = (all_loop_hg38_3loops[0:20])['loop_length'].map(color_dict) # ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
g = sns.barplot(data = all_loop_hg38_3loops[0:20], x = 'loop', y = 'value', palette = palette)
#g=sns.pointplot(data = all_loop_hg38_4333[0:20], x='loop', y='value', hue='loop_length')
g.set_xlabel("Loops", fontsize=20)
g.set_ylabel("Counts of each specific loop", fontsize = 20)
g.tick_params(labelsize=16)
g.axes.set_title("The distributiion of the frist 3 loops for all G4 ", fontsize=24)
handles, labels = g.get_legend_handles_labels()
g.legend()
plt.tight_layout()
plt.savefig('/uufs/chpc.utah.edu/common/home/u0674686/Desktop/hg38_3loops_distribution', dpi = 500)

hg38_loops = hg38[['g_runs', 'loop1', 'loop2', 'loop3', 'loop4', 'loop5']]
melt_hg38_loops = pd.melt(hg38_loops, id_vars = ['g_runs'], var_name = 'all_loop') ## melt all three loops to one column
all_loop_hg38_loops=(melt_hg38_loops.value.value_counts()).to_frame()
all_loop_hg38_loops['loop']=all_loop_hg38_loops.index
all_loop_hg38_loops['loop_length'] = [len(x) for x in all_loop_hg38_loops.loop]


all_loop_hg38_loops=all_loop_hg38_loops[all_loop_hg38_loops.loop!='None']

unique(all_loop_hg38_loops.loop_length[0:100])

-

melt_hg38_loops[melt_hg38_loops.value=='AGA']#

hg38[hg38.g_runs=='4:4:6:6:3:3:'][hg38.loop1=='ACTGTTGT']

hg38_3333[hg38_3333.loop1=='A'][hg38_3333.loop2=='A'].shape## 1480 same G4 with single A loop

hg38_3333[hg38_3333.loop1=='T'][hg38_3333.loop2=='T'][hg38_3333.loop3=='T'] # 199 with single T loop

hg38_3333[hg38_3333.loop1=='C'][hg38_3333.loop2=='C'][hg38_3333.loop3=='C'].shape ## 22 with single C loop

hg38.columns

hg38=hg38.drop('Unnamed: 0', axis=1) ## drop the unwanted column


hg38_count=Counter(hg38.to_records(index=False).tolist())

hg38_count.most_common()

hg38.drop_duplicates()


