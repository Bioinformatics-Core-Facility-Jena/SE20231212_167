from sklearn import decomposition
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

workingDirectory = '' #CHANGE THIS VARIABLE TO THE DESIRED WORKING DIRECTORY

##################################################
##### PCA based on the count data of all samples
##################################################

sampleMetaData = pd.read_csv(f'{workingDirectory}/sample_information.csv', sep=';')

countDF = pd.read_csv(f'{workingDirectory}/counts/all_samples_normalized_counts.csv', index_col=0)

#################################################################
##### PCA based on all Genes
countDFst =  StandardScaler().fit_transform(countDF) 
pcaOut = decomposition.PCA().fit(countDFst)

loadings_df = pd.DataFrame.from_dict(dict(zip(["PC"+str(i) for i in list(range(1, pcaOut.n_features_in_+1))], pcaOut.components_)))
loadings_df['treatment'] = sampleMetaData['Group']
loadings_df['ID'] = sampleMetaData['Sample_ID']

colors = sns.color_palette('Set1')
colors.extend(sns.color_palette('Pastel1'))


fig, ax = plt.subplots(1,2, figsize=(20, 12))
sns.scatterplot(data=loadings_df, x="PC1", y="PC2", hue='treatment', palette=colors, s=100, alpha=0.8, legend='full', edgecolor='black', ax=ax[0])
ax[0].set_xlabel(f"PC1\n({round(pcaOut.explained_variance_ratio_[0]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
ax[0].set_ylabel(f"PC2\n({round(pcaOut.explained_variance_ratio_[1]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
sns.scatterplot(data=loadings_df, x="PC1", y="PC2", s=250, alpha=0.8, legend='full', edgecolor='black', ax=ax[1])
ax[1].set_xlabel(f"PC1\n({round(pcaOut.explained_variance_ratio_[0]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
ax[1].set_ylabel(f"PC2\n({round(pcaOut.explained_variance_ratio_[1]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
d ={x:[float(loadings_df[loadings_df['ID'] == x]['PC1']), float(loadings_df[loadings_df['ID'] == x]['PC2'])] for x in loadings_df['ID']}
for idx,key in enumerate(d):
    ax[1].text(d[key][0], d[key][1], str(idx+1), ha='center', va='center', fontsize=8, weight='bold', color='white')
    ax[1].plot(d[key][0], d[key][1], label=f'{idx+1} -- {key}', alpha=0.0)
plt.legend(loc = 2, bbox_to_anchor = (1,1), title='Samples')
fig.get_figure()
plt.suptitle('PCA based on the expression values of all genes', fontsize=30)
fig.tight_layout()
fig.savefig(f'{workingDirectory}/pca_all_samples.pdf', dpi=300)
fig.clear()


#################################################################
##### PCA based on the top 500 variant genes
countDFst = pd.DataFrame(countDFst, columns = countDF.columns.values)
countDFst.set_index(countDF.index, inplace = True)

indexVarPairs = []
for rowname in countDFst.index:
    indexVarPairs.append((rowname, np.var(countDFst.loc[rowname])))

pcaOut = decomposition.PCA().fit(countDFst.loc[[x[0] for x in sorted(indexVarPairs, key=lambda x: x[1], reverse=True)[:500]]])
loadings_df = pd.DataFrame.from_dict(dict(zip(["PC"+str(i) for i in list(range(1, pcaOut.n_features_+1))], pcaOut.components_)))
loadings_df['treatment'] = sampleMetaData['Group']
loadings_df['ID'] = sampleMetaData['Sample_ID']

colors = sns.color_palette('Set1')
colors.extend(sns.color_palette('Pastel1'))

fig, ax = plt.subplots(1,2, figsize=(20, 12))
sns.scatterplot(data=loadings_df, x="PC1", y="PC2", hue='treatment', palette=colors, s=100, alpha=0.8, legend='full', edgecolor='black', ax=ax[0]) 
ax[0].set_xlabel(f"PC1\n({round(pcaOut.explained_variance_ratio_[0]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
ax[0].set_ylabel(f"PC2\n({round(pcaOut.explained_variance_ratio_[1]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
sns.scatterplot(data=loadings_df, x="PC1", y="PC2", s=250, alpha=0.8, legend='full', edgecolor='black', ax=ax[1])
ax[1].set_xlabel(f"PC1\n({round(pcaOut.explained_variance_ratio_[0]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
ax[1].set_ylabel(f"PC2\n({round(pcaOut.explained_variance_ratio_[1]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
d ={x:[float(loadings_df[loadings_df['ID'] == x]['PC1']), float(loadings_df[loadings_df['ID'] == x]['PC2'])] for x in loadings_df['ID']}
for idx,key in enumerate(d):
    ax[1].text(d[key][0], d[key][1], str(idx+1), ha='center', va='center', fontsize=8, weight='bold', color='white')
    ax[1].plot(d[key][0], d[key][1], label=f'{idx+1} -- {key}', alpha=0.0)
plt.legend(loc = 2, bbox_to_anchor = (1,1), title='Samples')
fig.get_figure()
plt.suptitle('PCA based on the top 500 genes with highest expression variance', fontsize=30)
fig.tight_layout()
fig.savefig(f'{workingDirectory}/pca_all_samples_top500.pdf', dpi=300)
fig.clear()

#################################################################
##### PCA based on the top 100 variant genes
countDFst = pd.DataFrame(countDFst, columns = countDF.columns.values)
countDFst.set_index(countDF.index, inplace = True)

indexVarPairs = []
for rowname in countDFst.index:
    indexVarPairs.append((rowname, np.var(countDFst.loc[rowname])))

pcaOut = decomposition.PCA().fit(countDFst.loc[[x[0] for x in sorted(indexVarPairs, key=lambda x: x[1], reverse=True)[:100]]])
loadings_df = pd.DataFrame.from_dict(dict(zip(["PC"+str(i) for i in list(range(1, pcaOut.n_features_+1))], pcaOut.components_)))
loadings_df['treatment'] = sampleMetaData['Group']
loadings_df['ID'] = sampleMetaData['Sample_ID']

colors = sns.color_palette('Set1')
colors.extend(sns.color_palette('Pastel1'))

fig, ax = plt.subplots(1,2, figsize=(20, 12))
sns.scatterplot(data=loadings_df, x="PC1", y="PC2", hue='treatment', palette=colors, s=100, alpha=0.8, legend='full', edgecolor='black', ax=ax[0])
ax[0].set_xlabel(f"PC1\n({round(pcaOut.explained_variance_ratio_[0]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
ax[0].set_ylabel(f"PC2\n({round(pcaOut.explained_variance_ratio_[1]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
sns.scatterplot(data=loadings_df, x="PC1", y="PC2", s=250, alpha=0.8, legend='full', edgecolor='black', ax=ax[1])
ax[1].set_xlabel(f"PC1\n({round(pcaOut.explained_variance_ratio_[0]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
ax[1].set_ylabel(f"PC2\n({round(pcaOut.explained_variance_ratio_[1]*100, 2)}% explained Variance)", fontdict={'weight':'bold', 'size':12}) 
d ={x:[float(loadings_df[loadings_df['ID'] == x]['PC1']), float(loadings_df[loadings_df['ID'] == x]['PC2'])] for x in loadings_df['ID']}
for idx,key in enumerate(d):
    ax[1].text(d[key][0], d[key][1], str(idx+1), ha='center', va='center', fontsize=8, weight='bold', color='white')
    ax[1].plot(d[key][0], d[key][1], label=f'{idx+1} -- {key}', alpha=0.0)
plt.legend(loc = 2, bbox_to_anchor = (1,1), title='Samples')
fig.get_figure()
plt.suptitle('PCA based on the top 100 genes with highest expression variance', fontsize=30)
fig.tight_layout()
fig.savefig(f'{workingDirectory}/pca_all_samples_top100.pdf', dpi=300)
fig.clear()


##### PCA Grid
def pca_grid(n=4):
    fig, ax = plt.subplots(n-1, n-1, figsize=(n*2.25, n*1.5))
    for i,j in [(i,j+1) for i in range(n-1) for j in range(i, n-1)]:
        if i == 0 and j == 1:  # Create a legend for the first subplot only
            pcaPlot = sns.scatterplot(data=loadings_df, x=f"PC{j+1}", y=f"PC{i+1}", hue='treatment', palette=colors, s=n*4, alpha=0.8, legend='brief', edgecolor='black', ax=ax[i][j-1])
        else:
            pcaPlot = sns.scatterplot(data=loadings_df, x=f"PC{j+1}", y=f"PC{i+1}", hue='treatment', palette=colors, s=n*4, alpha=0.8, legend=False, edgecolor='black', ax=ax[i][j-1])
        ax[i][j-1].set_ylabel('')
        ax[i][j-1].set_xlabel('')
        if i == 0:
            ax[i][j-1].set_title(f'PCA {j+1}', fontdict={'weight':'bold', 'size':12})
        if i == j-1:
            ax[i][j-1].set_ylabel(f'PCA {i+1}', fontdict={'weight':'bold', 'size':12}, rotation="horizontal", labelpad=30)
    for i in range(n-1):
        for j in range(n-1):
            if i > j:
                ax[i][j].remove()
                continue
    fig.tight_layout()
    handles, labels = ax[0][0].get_legend_handles_labels()
    legend = ax[0][0].legend(loc='lower left')
    legend.remove()
    fig.legend(handles, labels, loc='lower left')
    fig.savefig(f'{workingDirectory}/pca_grid.pdf', dpi=300)
    fig.clear()

pca_grid(5)
