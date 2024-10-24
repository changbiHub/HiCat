# %%
import os
import subprocess
import platform

# Get R_HOME
if platform.system() == 'Windows':
    r_home = r"C:/PROGRA~1/R/R-44~1.1" # Need to update this path for your system if you are using Windows
else:
    r_home = subprocess.check_output(['R', 'RHOME'], universal_newlines=True).strip()


# Get R binary path
r_path = os.path.join(r_home, 'bin')

# Update PATH
os.environ['R_HOME'] = r_home
os.environ['PATH'] = f"{r_path}:{os.environ['PATH']}"

print(f"R_HOME is set to: {os.environ['R_HOME']}")
print(f"PATH has been updated to include: {r_path}")

# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import seaborn as sns

# %%
import pyreadr
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# %%
from catboost import Pool,CatBoostClassifier
import umap
from sklearn.cluster import DBSCAN
from sklearn.metrics import accuracy_score,precision_score, recall_score, f1_score, confusion_matrix
import argparse

# %%
def plot_confusion_matrix(y_true, y_pred):
    # Create confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    
    # Get unique labels
    labels = np.unique(np.concatenate((y_true, y_pred)))
    
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot heatmap
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=labels, yticklabels=labels, ax=ax)
    
    # Set labels and title
    ax.set_xlabel('Predicted')
    ax.set_ylabel('True')

    
    # Display the plot
    plt.tight_layout()


# %%
# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Create an argument parser
parser = argparse.ArgumentParser(description='Script to perform cell annotation.')

# Add an argument for the base path
parser.add_argument('--base_path', type=str, default=os.path.join(script_dir,"..",'data'), help='Base path for train and test data')
parser.add_argument('--outputFolderName', default=os.path.join(script_dir,"..",'outputs'), type=str, help='folder name for output')

# Parse the command line arguments
args = parser.parse_args()

# Get the base path from the command line arguments
base_path = args.base_path

# Get the folder name from the base path
folder_name = os.path.basename(os.path.normpath(base_path))
output_folder_name = args.outputFolderName

print(f"Script directory: {script_dir}")
print(f"Base path: {base_path}")


train_path = os.path.join(base_path, "train.rds")
test_path = os.path.join(base_path, "test.rds")
output_path = os.path.abspath(os.path.join(script_dir, "..", output_folder_name))

if not os.path.exists(output_path):
    os.makedirs(output_path)

# %%
# Activate automatic conversion between R and pandas DataFrames
pandas2ri.activate()

# Import required R libraries
base = importr('base')
seurat = importr('Seurat')
harmony = importr('harmony')
dplyr = importr('dplyr')
ggplot2 = importr('ggplot2')
cowplot = importr('cowplot')

# Read the R script
with open(os.path.join(script_dir, 'harmony.r'), 'r') as file:
    r_script = file.read()

# Replace the placeholders in the R script with actual paths
r_script = r_script.replace('train_path', f'"{train_path}"'.replace("\\", "\\\\"))
r_script = r_script.replace('test_path', f'"{test_path}"'.replace("\\", "\\\\"))
r_script = r_script.replace('output_path', f'"{output_path}"'.replace("\\", "\\\\"))

# Print the modified R script for debugging
print(r_script)

# Execute the R script
try:
    harmony_result = robjects.r(r_script)
except Exception as e:
    print(f"Error executing R script: {e}")
    raise

# Convert R data frames to pandas DataFrames
train_df = pandas2ri.rpy2py(harmony_result[0])
test_df = pandas2ri.rpy2py(harmony_result[1])

train_df["cellType"] = pyreadr.read_r(train_path)[None]["yy"].values
train_df["train"] = True

test_df["cellType"] = pyreadr.read_r(test_path)[None]["yy"].values
test_df["train"] = False

df = pd.concat([train_df, test_df])

# %%
allTypes = df.cellType.unique()
knownTypes = df[df.train].cellType.unique()
unseenTypes = list(set(allTypes)-set(knownTypes))

# %%
# Create an instance of UMAP
# umap_model = umap(n_components=2)
umap_model = umap.UMAP(n_components=2)

# Transform count_pca using UMAP
umap_embedding = umap_model.fit_transform(df.drop(["cellType","train"],axis=1))

# Create a DataFrame with UMAP embedding and cluster labels
umap_df = pd.DataFrame({'UMAP_1': umap_embedding[:, 0], 'UMAP_2': umap_embedding[:, 1], 'cell type': df.cellType.values})
# Plot the UMAP embedding with Plotly
fig = px.scatter(umap_df, x='UMAP_1', y='UMAP_2', color='cell type', width=800, height=500)
fig.write_html(os.path.join(output_path,"umap_cellTypes.html"))

# %%
# Create an instance of DBSCAN
# dbscan_model = DBSCAN(eps=0.5,min_samples=20)
dbscan_model = DBSCAN()

# Apply DBSCAN on umap_embedding
dbscan_labels = dbscan_model.fit_predict(umap_embedding)

# Create column names dynamically based on the number of UMAP dimensions
umap_columns = [f'UMAP_{i+1}' for i in range(umap_embedding.shape[1])]

# Create a dictionary with UMAP dimensions
umap_dict = {col: umap_embedding[:, i] for i, col in enumerate(umap_columns)}

# Add the DBSCAN labels to the dictionary
umap_dict['cluster'] = dbscan_labels.astype(str)

# Create the DataFrame
umap_df = pd.DataFrame(umap_dict)

# Plot the UMAP embedding with Plotly
fig = px.scatter(umap_df, x='UMAP_1', y='UMAP_2', color='cluster', width=800, height=500)

fig.write_html(os.path.join(output_path,"umap_clusters.html"))

# %%
# Convert umap_embedding and dbscan_labels to DataFrames
pca_df = df.drop(["cellType","train"],axis=1).reset_index(drop=True)

# Concatenate the DataFrames
concatenated_df = pd.concat([umap_df, pca_df], axis=1)

X_train = concatenated_df[df.train.values]
X_test = concatenated_df[~df.train.values]
y_train = df[df.train.values].cellType
y_test = df[~df.train.values].cellType

# %%
# Create an instance of CatBoostClassifier
catboost_model = CatBoostClassifier(random_seed=42)

# Fit the model on the filtered samples with new features
catboost_model.fit(X_train.values, y_train.values.flatten())

# %%
# Get the predicted labels for the test set
y_pred = catboost_model.predict(X_test.values).flatten()
y_pred_proba = catboost_model.predict_proba(X_test.values)
confidence = y_pred_proba.max(axis=1)

# %%
y_test_mask = y_test.copy()
y_test_mask[y_test.isin(unseenTypes)] = "others"

# %%
# Assuming concatenated_df is your DataFrame
# and cat_features is a list of categorical feature indices or names

pool = Pool(data=concatenated_df.values)
feature_importances = catboost_model.get_feature_importance(data=pool)
# Normalize feature importances
normalized_importances = feature_importances / np.sum(feature_importances)
pd.DataFrame(normalized_importances, index=concatenated_df.columns, columns=['importance']).sort_values(by='importance', ascending=False).to_csv(os.path.join(output_path,"feature_importance.csv"))

# %%
rank_confidence = pd.DataFrame({"cluster": X_test["cluster"].values, "confidence": confidence}).groupby("cluster").mean().sort_values(by="confidence",ascending=True)
num_confident_cluster = len(rank_confidence)-np.nanargmax(rank_confidence.diff())
rank_confidence_th_idx = len(knownTypes)
unconfident_cluster = rank_confidence.iloc[:-rank_confidence_th_idx].index.values
confidence_th = rank_confidence.iloc[np.nanargmax(rank_confidence.diff())].values[0]
rank_confidence.reset_index(inplace=True)
rank_confidence["threshold"] = False
rank_confidence.loc[rank_confidence["confidence"]==confidence_th,"threshold"] = True
rank_confidence.to_csv(os.path.join(output_path,"rank_confidence.csv"),index=False)

# %%
y_pred_final = y_pred.copy()
y_pred_final[confidence<confidence_th] = X_test[confidence<confidence_th]["cluster"]

# %%
knownTypes_pred = list(set(y_pred[np.isin(y_pred,knownTypes)&(confidence>=confidence_th)]))
knownTypes_pred.sort()
knownTypes.sort()
unseenTypes.sort()
index_order = list(knownTypes)+list(unseenTypes)

# %%
result_final = pd.DataFrame({"y_true": y_test.values, "y_pred": y_pred_final})
result_final["count"] = 1
result_final = result_final.groupby(["y_true", "y_pred"]).count().reset_index()
pivot_table = pd.pivot_table(result_final, values='count', index='y_true', columns='y_pred', fill_value=0)
# Create a categorical index with the desired order
pivot_table.index = pd.Categorical(pivot_table.index, categories=index_order, ordered=True)

# Sort the pivot table based on the categorical index
sorted_pivot_table = pivot_table.sort_index()
sorted_pivot_table[list(knownTypes_pred)+list(set(sorted_pivot_table.columns)-set(knownTypes))].to_csv(os.path.join(output_path,"confusion_matrix.csv"))

# %%
y_pred_final_masked = y_pred_final.copy()
y_pred_final_masked[~np.isin(y_pred_final, knownTypes)] = "others"
plot_confusion_matrix(y_test_mask, y_pred_final_masked)
plt.savefig(os.path.join(output_path,"confusion_matrix.png"))

# %%
pd.DataFrame({"accuracy":accuracy_score(y_test_mask, y_pred_final_masked), 
              "high_confidence_accuracy":accuracy_score(y_test_mask[confidence>=confidence_th], y_pred_final_masked[confidence>=confidence_th]), 
              "precision_unseen":precision_score(y_test_mask, y_pred_final_masked, average=None, labels=["others"])[0], 
              "recall_unseen":recall_score(y_test_mask, y_pred_final_masked, average=None, labels=["others"])[0], 
              "f1_unseen":f1_score(y_test_mask, y_pred_final_masked, average=None, labels=["others"])[0]},index=[0]).to_csv(os.path.join(output_path,"metrics.csv"),index=False)

pd.DataFrame({"y_true":y_test.values,"y_pred_supervised":y_pred,"y_pred_combined":y_pred_final}).to_csv(os.path.join(output_path,"y_pred.csv"),index=False)