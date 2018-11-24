from __future__ import division
import numpy as np
import random
import re
import pandas as pd
import matplotlib.pyplot as plt
#%%
def random_seq(gl=1000000, GC=50):
    """a function to generate a random sequences that the length are gl, GC percentage are GC"""
    base=['A', 'T', 'C', 'G']
    base_distr=[(1-float(GC/100))/2, (1-float(GC/100))/2, float(GC/100)/2, float(GC/100)/2]
    seq=''.join(np.random.choice(base,gl,p=base_distr))
    return seq
random_seq().count('A')
r1=random_seq(gl=150, GC=60)
print(r1, (r1.count('G')+ r1.count('C'))/len(r1))


def n_random_seq(n=10, GCs=list(range(20,90,10)), gl=1000000):
    """genereate n random sequences for each GC percentage in GCs, return a dictionary"""
    all_seq={}
    for i in GCs:
        all_seq[i]=[random_seq(gl=gl,GC=i) for x in range(n)]
    return all_seq



def g4_counter(seq):
    """take a sequences count the number of G4
    The default n is 3, which is looking for more than 4 tracts of 3 Gs
    when n = 4, it is looking for PQS with a spare tire"""
    g4=r'((G){3,}\w{1,12}){3,}(G){3,}'
    c4=r'((C){3,}\w{1,12}){3,}(C){3,}'
    n_g4=len(tuple(re.finditer(g4,seq)))
    n_c4=len(tuple(re.finditer(c4,seq)))
    g5=r'((G){3,}\w{1,12}){4,}(G){3,}'
    c5=r'((C){3,}\w{1,12}){4,}(C){3,}'
    n_g5=len(tuple(re.finditer(g5,seq)))
    n_c5=len(tuple(re.finditer(c5,seq)))
    all_g4=n_g4+n_c4
    all_g5=n_g5 + n_c5
    try:
        spare_percent =float(all_g5)/float(all_g4) 
    except ZeroDivisionError:
        spare_percent = 0
    return all_g4, spare_percent

# In[79]:

#np.mean([len(tuple(re.finditer(g4,n_random_seq()[70][3]))) for i in range(100)]) ## count by finditer


# In[95]:

def g4_counter_from_dict(g4_dict):
    """take a g4 dictionary and count the number of G4s for each GC percentage, return a dictionary
    the default is looking for 4 tracks of 3 Gs or more (n = 3)"""
    g4_n={}
    g4_spare = {}
    for k,v in g4_dict.items():
        tt=[]
        spare=[]
        for j in v:
            g4_number, spare_percent = g4_counter(j)
            tt.append(g4_number)
            spare.append(spare_percent)
        g4_n[k]=tt
        g4_spare[k] = spare
    return g4_n, g4_spare


"""g4_5=n_random_seq(n=10, GCs=list(range(20,80,5)))
tt=g4_counter_from_dict(g4_5)

pd_tt=pd.DataFrame.from_dict(tt) ## save data to pandas


# In[106]:

pd_tt.to_csv("/home/jovyan/work/random_g4_at_different_GC_content.csv")## save to csv file"""


g4_10=n_random_seq(n=50, GCs=list(range(20,80,10)), gl=2000000) ## every 2 percent of GC, calculate random G4s
tt_10_gl_2m, tt_10_spare=g4_counter_from_dict(g4_10)## a dictionary
#tt_10_gl_2m_spare = g4_counter_from_dict(g4_10, n=4)
pd_tt_all = pd.DataFrame.from_dict(tt_10_gl_2m)
pd_tt_spare_percent = pd.DataFrame.from_dict(tt_10_spare)
print(pd_tt_spare_percent)
pd_tt_all.to_csv("~/Documents/bacteria_G4/G4_count_random_2m_20_80_all.csv")
pd_tt_spare_percent.to_csv("~/Documents/bacteria_G4/G4_count_random_2m_20_80_spare_percent.csv")
pd_tt_spare_percent.boxplot()
plt.show()
pd_tt_all.boxplot()
plt.show()

tt_10_gl_2m_density={}
tt_10_gl_2m_spare_density = {}
for k, v in tt_10_gl_2m.items():
    tt_10_gl_2m_density[k]=[i/2 for i in v] ## calculate G4 density per 1 Mb
for k, v in tt_10_gl_2m_spare.items():
    tt_10_gl_2m_spare_density[k]=[i/2 for i in v] ## calculate G4 density per 1 Mb

pd_2m_2=pd.DataFrame.from_dict(tt_10_gl_2m_density)
pd_2m_2_spare=pd.DataFrame.from_dict(tt_10_gl_2m_spare_density)

pd_2m_2.plot.box()

pd_2m_2.to_csv("~/Documents/bacteria_G4/G4_count_random_2m_20_80_all.csv")
pd_2m_2_spare.to_csv("~/Documents/bacteria_G4/G4_count_random_2m_20_80_spare.csv")