
#%%
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import pandas as pd
import numpy as np
import pickle as pkl
from scipy.io import arff

#%% [markdown]
# # Data Exploration
# 
# We are analyzing a Breast Cancer Data Set with 286 instances of real patient data from the Institute of Oncology University Medical Center in Ljubljana, Yugoslavia. The data set is available at the UCI repository
# (http://archive.ics.uci.edu/ml/datasets/Breast+Cancer) and contains the following variables:
# 
# Class (recurrence): This is our label, recurrence-events vs no-recurrence-events
# 1. Age: age (in years at last birthday) of the patient at the time of diagnosis;
# 2. Menopause: whether the patient is pre- or postmenopausal at time of diagnosis;
# 3. Tumor_Size: the greatest diameter (in mm) of the excised tumor;
# 4. Inv-nodes: the number (range 0 - 39) of axillary lymph nodes that contain metastatic breast cancer visible on histological examination;
# 5. Node-caps: if the cancer does metastasise to a lymph node, although outside the original site of the tumor it may remain “contained” by the capsule of the lymph node. However, over time, and with a more aggressive disease, the tumor may replace the lymph node and then penetrate the capsule, allowing it to invade the surrounding tissues;
# 6. Deg-malig (Degree of malignancy): Histological Grade (range 1-3) of the tumor (degree of abnormality displayed by the cells). Tumors that are grade 1 predominantly consist of cells that, while neoplastic, retain many of their usual characteristics.
# Grade 3 tumors predominately consist of cells that are highly abnormal;
# 7. Breast: breast cancer may obviously occur in either
# breast;
# 8. Breast-quad (Breast quadrant): the breast may be divided into four quadrants, using the nipple as a central point;
# 9. Irradiat (Irradiation): radiation therapy is a treatment that uses high-energy x-rays to destroy cancer cells. 

#%%
#Data read in and processing:
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer/breast-cancer.data"
names = ['recurrence', 'Age', 'Menopause', 'Tumor_Size', 'Inv-nodes', 'Node-caps', 'Deg-malig', 'Breast', 'Breast_quad', 'Irradiat']
bc = pd.read_csv(url, names=names)
print(bc.shape)
bc.head()
#array = dataframe.values


#%%
#Initial variant and data exploration
for i in bc.columns:
    print(i, bc[i].value_counts())
# find column 5 and 8 has '?', which is the missing value


#%%
# remove 9 instances with missing values:
bc=bc[bc['Node-caps']!='?'][bc['Breast_quad']!='?']
bc.shape


#%%
#Initial variant and data exploration
print("\n \t The data frame has %d rows and %d columns. \n"%(bc.shape[0], bc.shape[1]))
bc.info()

#%% [markdown]
# From the information above, we can see that except for "deg-malig", all other features are of type object and have 0 non-null numbers. Because many algorithms cannot work with categorical data, we need to transform our features of type "object" and encode the data. We are going to use One Hot encoding for features 7 and 8. And we are using Integer encoding for label variable and our features 1-5,9. We are not encoding "deg-malig" because it is already numeric. 
# 
# .

#%%
def int_encoder(intgr_encoding): 
    # integer encode
    label_encoder = LabelEncoder()
    integer_encoded = label_encoder.fit_transform(intgr_encoding)
    return integer_encoded


#%%
# integer encode
bc['recurrence'] = int_encoder(intgr_encoding = bc['recurrence'])
bc['Age'] = int_encoder(intgr_encoding =  bc['Age'])
bc['Menopause'] = int_encoder(intgr_encoding = bc['Menopause'])
bc['Tumor_Size'] = int_encoder(intgr_encoding = bc['Tumor_Size'])
bc['Inv-nodes'] = int_encoder(intgr_encoding = bc['Inv-nodes'])
bc['Node-caps'] = int_encoder(intgr_encoding = bc['Node-caps'])
bc['Irradiat'] = int_encoder(intgr_encoding = bc['Irradiat'])
bc.head()


#%%
# one hot encoding
bc = pd.get_dummies(bc, columns=['Breast','Breast_quad'])
bc.head()


#%%
bc['recurrence'] = bc['recurrence'].astype('category')


#%%
print("\n \t The data frame has %d rows and %d columns. \n"%(bc.shape[0], bc.shape[1]))
bc.dtypes


#%%



#%%



