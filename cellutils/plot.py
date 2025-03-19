import os
import re

import numpy as np
import pandas  as pd

import matplotlib.pyplot as plt
import seaborn as sns

from synergy.combination.bliss import Bliss
from synergy.utils.plots import aggregate_replicates, is_on_grid

def plot_plate(df:pd.DataFrame, 
               cond="COND", index='index', score='score', count=None, wellid='Image_Metadat_WellID',
               cutoff=None, data_cols=[], title='well plot',save_img=False,img_name='my_plot.png', 
               save_hits=False, hits_name='my_hits.csv', cmpd='CMPD'):
    """
    Plot a 384 well plate's scores, coloring on compounds and contorls
    Args:
        df (pd.DataFrame): well level dataframe 
        cond (str, optional): Condition of well (PC, NC, CMPD) Defaults to "COND".
        index (str, optional): Well number (1-384). Defaults to 'index'.
        score (str, optional): Y-axis Defaults to 'score'.
        count (_type_, optional): size feature. Defaults to None.
        wellid (str, optional): WellId used for finding hits Defaults to 'Image_Metadat_WellID'.
        cutoff (_type_, optional): set cutoff hits mean+(3*std). Defaults to None.
        data_cols (list, optional): additional columns hit csv Defaults to [].
        title (str, optional): title of plot Defaults to 'well plot'.
        save_img (bool, optional): save figure Defaults to False.
        img_name (str, optional): figure name to save to Defaults to 'my_plot.png'.
        save_hits (bool, optional): save hits csv Defaults to False.
        hits_name (str, optional): name of hits file Defaults to 'my_hits.csv'.
    """
    cols = [index,wellid]
    std = df.loc[df[cond]==cmpd, score].std()
    mean = df.loc[df[cond]==cmpd, score].mean()
    if cutoff is None:
        cutoff = mean + (3*std)
    pnts = df.loc[(df[score]>=cutoff)&(df[cond]==cmpd), data_cols+cols]

    ax = sns.scatterplot(data=df, x=index, y=score, size=count, hue=cond)
    plt.axhline(y=cutoff, c='r')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1,1))

    for i, point in pnts.iterrows():
        ax.text(point[index]+2, point[score]+0.01, str(point[wellid]))

    plt.title(title)

    if save_img:
        plt.savefig(img_name, format='png', bbox_inches="tight")
    print(pnts[wellid].unique())
    if save_hits:
        pnts.to_csv(hits_name, index=False)

    return pnts

def plot_reg_facet(df:pd.DataFrame, x='concentration', y='count', 
                   x_size=4, y_size=4, cmpd_col='Compound', 
                   plot_name='my_cool_plot',xlabel='Log Concentration', ylabel='Percent Viabilty'):
    """
    plot_reg_facet
    Plots a facet grid of regression plots
    """ 
    cmpds = df[cmpd_col].unique().tolist()
    fig, ax = plt.subplots(x_size, y_size, sharex=True, sharey=True)
    fig.set_figheight(30)
    fig.set_figwidth(30)
    
    idx = 0
    
    for row in range(0, x_size):
        for col in range(0, y_size):
            if idx > len(cmpds)-1:
                break
            c = cmpds[idx]
            sns.regplot(ax=ax[row, col], x=x, y=y, data=df.loc[df[cmpd_col]==c])
            ax[row,col].set_title(c, fontsize=15, fontweight="bold")
            ax[row,col].set_xlabel(xlabel)
            ax[row,col].set_ylabel(ylabel)
            idx+=1
    
    plt.savefig(f"{plot_name}_regression.svg", format='SVG', dpi=2000)
    
def synergy_plot(df, drugs, cell_line, agaisnt, bliss=False, output_dir=''):
    for dr in drugs:
        print(cell_line, dr)
        a = df.loc[df['Drug1']==dr]
        title = f'{agaisnt} vs {dr}'
        xlabel = agaisnt
        d1 = np.array(a['Conc1'].values, copy=True, dtype=np.float64)
        d2 = np.array(a['Conc2'].values, copy=True, dtype=np.float64)
        vals = np.array(a['score'].values, copy=True, dtype=np.float64)
        
        if bliss:
            model = Bliss()
            vals = model.fit(d1,d2, vals)
        
        sorted_indices = np.lexsort((d1,d2))
        D1 = d1[sorted_indices]
        D2 = d2[sorted_indices]
        vals = vals[sorted_indices]
        D_unique, vals = aggregate_replicates(np.vstack((D1,D2)).T, vals, aggfunc=np.median)
        if not is_on_grid(D_unique):
            raise ValueError("plot_surface_plotly() requires d1, d2 to represent a dose grid")
        D1 = D_unique[:, 0]
        D2 = D_unique[:, 1]
        dy = pd.DataFrame(data={
            'conc1':D_unique[:,0],
            'conc2':D_unique[:,1],
            'response':vals
        })
        df_p = dy.pivot(index='conc1', columns='conc2', values='response')
        ax = sns.heatmap(data=df_p, annot=True,fmt='.2f')
        ax.set(xlabel=xlabel,ylabel=dr, title=title)
        ax.invert_yaxis()
        if bliss:
            fname = os.path.join(output_dir, f"bliss_{cell_line}_{dr}.svg")
        else:
            fname = os.path.join(output_dir, f"{cell_line}_{dr}.svg")
        fig = ax.get_figure()
        plt.savefig(fname, format='SVG')
        plt.clf()
