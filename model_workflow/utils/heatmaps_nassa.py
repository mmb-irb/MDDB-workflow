import pathlib
import itertools

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import ListedColormap
from model_workflow.utils.nucleicacid import NucleicAcid

def get_axes(subunit_len, base):
    '''Create the correct labels according to the lenght of the subunit'''

    yaxis = [
        f"{base}{base}",
        f"{base}C",
        f"C{base}",
        "CC",
        f"G{base}",
        "GC",
        f"A{base}",
        "AC",
        f"{base}G",
        f"{base}A",
        "CG",
        "CA",
        "GG",
        "GA",
        "AG",
        "AA"]

    if subunit_len == 3:
        xaxis = f"G A C {base}".split()
        nucleotide_order = [b.join(f) for f in yaxis for b in xaxis]

    elif subunit_len == 4:
        xaxis = yaxis[::-1]
        nucleotide_order = [b.join(f) for f in yaxis for b in xaxis]

    elif subunit_len == 5:
        #xaxis = f"G A C {base}".split()
        xaxis = []
        yaxis = []
        nucl = ['A', 'C', 'G', base]
        nucleotide_order = [a + b + c + d + e for a in nucl for b in nucl for c in nucl for d in nucl for e in nucl]
        for i in nucleotide_order:
            xaxis.append(i[1:-1])
            yaxis.append(i[0] + "..." + i[-1:])
            #yaxis.append(i[:2] + "_" + i[-2:])
        yaxis = list(set(xaxis))

    elif subunit_len == 6:
        nucleotide_order= [f"{n1}{n2}{n3}" for n1 in yaxis for n2 in yaxis for n3 in yaxis]
        yaxis = []
        xaxis = []
        for i in nucleotide_order:
            yaxis.append(i[:2] + "_" + i[-2:])
            xaxis.append(i[2:4])
        yaxis = list(set(yaxis))
        xaxis = list(set(xaxis))

    return xaxis, yaxis, nucleotide_order


def reorder_labels_rotated_plot(df, subunit_name, tetramer_order):
    '''For the rotated plot: reorder the nucleotides labels matching them to their according values of the dataframe. '''
    
    if subunit_name == 'tetramer':
        sorted_index = dict(zip(tetramer_order, range(len(tetramer_order))))
        df["subunit_rank"] = df[subunit_name].map(sorted_index)
        df = df.sort_values(by="subunit_rank")
        df = df.drop("subunit_rank", axis=1)
        df = df.reset_index(drop=True)

    tetramer_order.sort()
    all_tetramers = pd.DataFrame({subunit_name: ["".join(tetramer) for tetramer in tetramer_order]})
    merged_df = pd.merge(all_tetramers, df, on=subunit_name, how='left')
    merged_df['subunit_rank'] = merged_df[subunit_name].map(lambda x: tetramer_order.index(x) if not pd.isna(x) else np.nan)
    merged_df = merged_df.sort_values(by='subunit_rank').reset_index(drop=True) 
    
    if subunit_name == 'hexamer':
        centre= len(merged_df[subunit_name][1])//2
        merged_df['yaxis'] = merged_df[subunit_name].apply(lambda x: x[:centre-2]+"...."+x[centre+2:])
        merged_df['xaxis'] = merged_df[subunit_name].apply(lambda x: x[centre-2:centre+2])
        merged_df = merged_df.sort_values(by=['yaxis','xaxis'], ascending=True)
        yaxiss = []
        xaxiss = []
        for i in range(len(merged_df)):
            xaxiss.append(merged_df[subunit_name][i][centre-2:centre+2])
            yaxiss.append(merged_df[subunit_name][i][:centre-2] + "...." + merged_df[subunit_name][i][centre+2:])

        xaxis1=[]
        for nucl in xaxiss:
            if nucl not in xaxis1:
                xaxis1.append(nucl)

        yaxis1 = []
        for nucl in yaxiss:
            if nucl not in yaxis1:
                yaxis1.append(nucl)
        df1 = merged_df

    if subunit_name == 'pentamer':
        centre= len(merged_df[subunit_name][1])%2
        merged_df['xaxis'] = merged_df[subunit_name].apply(lambda x: x[1:-1])
        merged_df['yaxis'] = merged_df[subunit_name].apply(lambda x: x[0]+"..."+x[-1:])
        merged_df = merged_df.sort_values(by=['yaxis','xaxis'], ascending=True)
        yaxiss = []
        xaxiss = []
        for i in range(len(merged_df)):
            xaxiss.append(merged_df[subunit_name][i][1:-1])
            yaxiss.append(merged_df[subunit_name][i][0] + "..." + merged_df[subunit_name][i][-1:])

        xaxis1=[]
        for nucl in xaxiss:
            if nucl not in xaxis1:
                xaxis1.append(nucl)

        yaxis1 = []
        for nucl in yaxiss:
            if nucl not in yaxis1:
                yaxis1.append(nucl)
        df1 = merged_df    

    # if subunit_name == 'pentamer':
    #     centre= len(merged_df[subunit_name][1])%2
    #     merged_df['xaxis'] = merged_df[subunit_name].apply(lambda x: x[:centre+1]+"_"+x[centre+2:])
    #     merged_df['yaxis'] = merged_df[subunit_name].apply(lambda x: x[centre+1])
    #     merged_df = merged_df.sort_values(by=['yaxis','xaxis'], ascending=True)
    #     yaxiss = []
    #     xaxiss = []
    #     for i in range(len(merged_df)):
    #         yaxiss.append(merged_df[subunit_name][i][centre+1])
    #         xaxiss.append(merged_df[subunit_name][i][:centre+1] + "_" + merged_df[subunit_name][i][centre+2:])

    #     xaxis1=[]
    #     for nucl in xaxiss:
    #         if nucl not in xaxis1:
    #             xaxis1.append(nucl)

    #     yaxis1 = []
    #     for nucl in yaxiss:
    #         if nucl not in yaxis1:
    #             yaxis1.append(nucl)
    #     df1 = merged_df   

    return df1, xaxis1, yaxis1

def reorder_labels_straight_plot(df, subunit_name, tetramer_order):
    '''For the straight plot: reorder the nucleotides labels matching them to their according values of the dataframe. '''

    if subunit_name == 'tetramer':
        sorted_index = dict(zip(tetramer_order, range(len(tetramer_order))))
        df["subunit_rank"] = df[subunit_name].map(sorted_index)
        df = df.sort_values(by="subunit_rank")
        df = df.drop("subunit_rank", axis=1)
        df = df.reset_index(drop=True)

    tetramer_order.sort()
    all_tetramers = pd.DataFrame({subunit_name: ["".join(tetramer) for tetramer in tetramer_order]})
    merged_df = pd.merge(all_tetramers, df, on=subunit_name, how='left')
    merged_df['subunit_rank'] = merged_df[subunit_name].map(lambda x: tetramer_order.index(x) if not pd.isna(x) else np.nan)
    merged_df = merged_df.sort_values(by='subunit_rank').reset_index(drop=True) 
    
    if subunit_name == 'hexamer':
        centre= len(merged_df[subunit_name][1])//2
        merged_df['yaxis'] = merged_df[subunit_name].apply(lambda x: x[:centre-1] + "_" + x[centre+1:])
        merged_df['xaxis'] = merged_df[subunit_name].apply(lambda x: x[centre-1:centre+1])
        merged_df = merged_df.sort_values(by=['yaxis','xaxis'], ascending=True)
        yaxiss = []
        xaxiss = []
        for i in range(len(merged_df)):
            xaxiss.append(merged_df[subunit_name][i][centre-1:centre+1])
            yaxiss.append(merged_df[subunit_name][i][:centre-1] + "_" + merged_df[subunit_name][i][centre+1:])

        xaxis=[]
        for nucl in xaxiss:
            if nucl not in xaxis:
                xaxis.append(nucl)

        yaxis = []
        for nucl in yaxiss:
            if nucl not in yaxis:
                yaxis.append(nucl)
        df = merged_df

    if subunit_name == 'pentamer':
        centre= len(merged_df[subunit_name][1])%2
        merged_df['yaxis'] = merged_df[subunit_name].apply(lambda x: x[1:-1])
        merged_df['xaxis'] = merged_df[subunit_name].apply(lambda x: x[0]+"..."+x[-1:])

        merged_df = merged_df.sort_values(by=['yaxis','xaxis'], ascending=True)
        yaxiss = []
        xaxiss = []
        for i in range(len(merged_df)):
            yaxiss.append(merged_df[subunit_name][i][1:-1])
            xaxiss.append(merged_df[subunit_name][i][0] + "..." + merged_df[subunit_name][i][-1:])

        xaxis=[]
        for nucl in xaxiss:
            if nucl not in xaxis:
                xaxis.append(nucl)

        yaxis = []
        for nucl in yaxiss:
            if nucl not in yaxis:
                yaxis.append(nucl)
        df = merged_df

    # if subunit_name == 'pentamer':
    #     centre= len(merged_df[subunit_name][1])%2
    #     merged_df['yaxis'] = merged_df[subunit_name].apply(lambda x: x[:centre+1] + "_" + x[centre+2:])
    #     merged_df['xaxis'] = merged_df[subunit_name].apply(lambda x: x[centre+1])
    #     merged_df = merged_df.sort_values(by=['yaxis','xaxis'], ascending=True)
    #     yaxiss = []
    #     xaxiss = []
    #     for i in range(len(merged_df)):
    #         xaxiss.append(merged_df[subunit_name][i][centre+1])
    #         yaxiss.append(merged_df[subunit_name][i][:centre+1] + "_" + merged_df[subunit_name][i][centre+2:])

    #     xaxis=[]
    #     for nucl in xaxiss:
    #         if nucl not in xaxis:
    #             xaxis.append(nucl)

    #     yaxis = []
    #     for nucl in yaxiss:
    #         if nucl not in yaxis:
    #             yaxis.append(nucl)
    #     df = merged_df

    return df, xaxis, yaxis

def arlequin_plot(
        df,
        global_mean,
        global_std,
        helpar,
        save_path,
        unit_name,
        unit_len,
        base,
        label_offset=0.5):
    
    xaxis, yaxis, tetramer_order = get_axes(unit_len, base)
    df, xaxis, yaxis = reorder_labels_straight_plot(df, unit_name, tetramer_order)
    df1, xaxis1, yaxis1 = reorder_labels_rotated_plot(df, unit_name, tetramer_order)

    sz1 = df["col1"].ravel()
    sz2 = df["col2"].ravel()
 
    if unit_name == 'hexamer':
        M = 4**2
        N = 4**(unit_len - 2)
        M_1 = 4 ** (unit_len - 2)
        N_1 = 4 ** 2
    
    if unit_name == 'pentamer':
        M = 4 * 4
        N = 4 ** 3
        M_1 = 4 ** 3
        N_1 = 4 * 4

    # STRAIGHT PLOT
    
    x = np.arange(M + 1)
    y = np.arange(N + 1)
    xs, ys = np.meshgrid(x, y)

    upper_triangle = [(i + j*(M+1), i+1 + j*(M+1), i+1 + (j+1)*(M+1))
                      for j in range(N) for i in range(M)]
    lower_triangle = [(i + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1))
                      for j in range(N) for i in range(M)]
    triang1 = Triangulation(xs.ravel(), ys.ravel(), upper_triangle)
    triang2 = Triangulation(xs.ravel(), ys.ravel(), lower_triangle)

    fig, axs = plt.subplots(
        1,      
        1,
        figsize=(8,18),
        dpi=300,
        tight_layout=True)

    colormap = plt.get_cmap("bwr", 3).reversed()
    colormap.set_bad(color="grey")
    img1 = axs.tripcolor(triang1, sz1, cmap=colormap, vmin=-1, vmax=1)
    _ = axs.tripcolor(triang2, sz2, cmap=colormap, vmin=-1, vmax=1)

    axs.grid()
    xlocs = np.arange(len(xaxis))
    ylocs = np.arange(len(yaxis))
    _ = axs.set_xticks(xlocs)
    _ = axs.set_xticklabels("")
    _ = axs.set_yticks(ylocs)
    _ = axs.set_yticklabels("")
    _ = axs.set_xticks(xlocs+label_offset, minor=True)
    _ = axs.set_xticklabels(xaxis, minor=True,fontsize=8)
    _ = axs.set_yticks(ylocs+label_offset, minor=True)
    _ = axs.set_yticklabels(yaxis, minor=True,fontsize=6)

    _ = axs.set_xlim(0, M)
    _ = axs.set_ylim(0, N)
    axs.set_title(helpar.upper())
    cbar = fig.colorbar(img1, ax=axs, ticks=[-1, 0, 1], shrink=0.4)
    cbar.ax.set_yticklabels([
        f"< {global_mean:.2f}-{global_std:.2f}",
        f"{global_mean:.2f}$\pm${global_std:.2f}",
        f"> {global_mean:.2f}+{global_std:.2f}"])

    file_path = pathlib.Path(save_path) / f"{helpar}.pdf"
    fig.savefig(fname=file_path, format="pdf")

    #  ROTATED PLOT

    sz1_1 = df1["col1"].ravel()
    sz2_1 = df1["col2"].ravel()

    x_1 = np.arange(M_1 + 1)
    y_1 = np.arange(N_1 + 1)
    xs_1, ys_1 = np.meshgrid(x_1, y_1)

    upper_triangle_1 = [(i + j*(M_1+1), i+1 + j*(M_1+1), i+1 + (j+1)*(M_1+1))
                      for j in range(N_1) for i in range(M_1)]
    lower_triangle_1 = [(i + j*(M_1+1), i+1 + (j+1)*(M_1+1), i + (j+1)*(M_1+1))
                      for j in range(N_1) for i in range(M_1)]
    triang1_1 = Triangulation(xs_1.ravel(), ys_1.ravel(), upper_triangle_1)
    triang2_1 = Triangulation(xs_1.ravel(), ys_1.ravel(), lower_triangle_1)

    fig_1, axs_1 = plt.subplots(
        1,      
        1,
        figsize=(22, 6), 
        dpi=300,
        tight_layout=True)

    colormap_1 = plt.get_cmap("bwr", 3).reversed()
    colormap_1.set_bad(color="grey")
    img1_1 = axs_1.tripcolor(triang1_1, sz1_1, cmap=colormap_1, vmin=-1, vmax=1)
    _ = axs_1.tripcolor(triang2_1, sz2_1, cmap=colormap_1, vmin=-1, vmax=1)

    axs_1.grid()
    xlocs_1 = np.arange(len(xaxis1))
    ylocs_1 = np.arange(len(yaxis1))
    _ = axs_1.set_xticks(xlocs_1)
    _ = axs_1.set_xticklabels("")
    _ = axs_1.set_yticks(ylocs_1)
    _ = axs_1.set_yticklabels("")
    _ = axs_1.set_xticks(xlocs_1+label_offset, minor=True)
    _ = axs_1.set_xticklabels(xaxis1, minor=True,fontsize=4, rotation=90)
    _ = axs_1.set_yticks(ylocs_1+label_offset, minor=True)
    _ = axs_1.set_yticklabels(yaxis1, minor=True,fontsize=6)
    _ = axs_1.set_xlim(0, M_1)
    _ = axs_1.set_ylim(0, N_1)
    axs_1.set_title(helpar.upper())
    cbar_1 = fig_1.colorbar(img1_1, ax=axs_1, ticks=[-1, 0, 1], shrink=0.4)
    cbar_1.ax.set_yticklabels([
        f"< {global_mean:.2f}-{global_std:.2f}",
        f"{global_mean:.2f}$\pm${global_std:.2f}",
        f"> {global_mean:.2f}+{global_std:.2f}"])

    file_path1 = pathlib.Path(save_path) / f"{helpar}_rotated.pdf"
    fig_1.savefig(fname=file_path1, format="pdf")
    return fig, axs, fig_1, axs_1 


def bconf_heatmap(df, fname, save_path, subunit_len, base="T", label_offset=0.05):
    print('subunit len: ', subunit_len)
    if subunit_len == 3:
        yaxis = [
        f"{base}{base}",
        f"{base}C",
        f"C{base}",
        "CC",
        f"G{base}",
        "GC",
        f"A{base}",
        "AC",
        f"{base}G",
        f"{base}A",
        "CG",
        "CA",
        "GG",
        "GA",
        "AG",
        "AA"]
        xaxis = f"G A C {base}".split()
        nucleotide_order = pd.DataFrame(
            [b.join(f) for f in yaxis for b in xaxis],
            columns=["trimer"])
        df = df.merge(nucleotide_order, how="right", on="trimer")
        fig, ax = plt.subplots()

    elif subunit_len == 4:
        xaxis = [
        "GG",
        "GA",
        "AG",
        "AA",
        "GC",
        f"G{base}",
        f"A{base}",
        "AC",
        "CA",
        f"{base}A",
        f"{base}G",
        "CG",
        "CC",
        f"C{base}",
        f"{base}C",
        f"{base}{base}"]
        yaxis = xaxis.copy()
        tetramer_order = pd.DataFrame(
            [b.join(f) for f in yaxis for b in xaxis],
            columns=["tetramer"])
        df = df.merge(tetramer_order, how="right", on="tetramer")
        fig, ax = plt.subplots()

    elif subunit_len == 5:
        raise ValueError('The length of the subunit 5 is not a valid option for bconf analyses. Try with 4 or 6.')

    elif subunit_len == 6:
        xaxis = [
        f"{base}{base}",
        f"{base}C",
        f"C{base}",
        "CC",
        f"G{base}",
        "GC",
        f"A{base}",
        "AC",
        f"{base}G",
        f"{base}A",
        "CG",
        "CA",
        "GG",
        "GA",
        "AG",
        "AA"]
        nucleotide_order= [f"{n1}{n2}{n3}" for n1 in xaxis for n2 in xaxis for n3 in xaxis]
        yaxis = []
        xaxis = []
        for i in nucleotide_order:
            xaxis.append(i[:2] + "_" + i[-2:])
            yaxis.append(i[2:4])
        yaxis = list(set(yaxis))
        xaxis = list(set(xaxis))
        # tetramer_order = pd.DataFrame(
        #     [b.join(f) for f in xaxis for b in yaxis],
        #     columns=["hexamer"])
        tetramer_order = pd.DataFrame(
            sorted(nucleotide_order),
            columns=["hexamer"])
        df = df.merge(tetramer_order, how="right", on="hexamer")
        fig, ax = plt.subplots(
            1,      
            1,
            figsize=(22, 12), 
            dpi=300,
            tight_layout=True)
    
    colormap = ListedColormap([
        "darkblue",
        "blue",
        "lightblue",
        "lightgreen",
        "lime",
        "orange",
        "red",
        "crimson"])
    colormap.set_bad(color="grey")

    # plot
    if subunit_len == 3:
        im = ax.imshow(df["pct"].to_numpy().reshape((16, 4)), cmap=colormap)
    if subunit_len == 4:
        im = ax.imshow(df["pct"].to_numpy().reshape((16, 16)), cmap=colormap)
    if subunit_len == 6:
        im = ax.imshow(df["pct"].to_numpy().reshape((16, 256)), cmap=colormap)
    plt.colorbar(im)
    # axes
    xlocs = np.arange(len(xaxis))
    ylocs = np.arange(len(yaxis))
    _ = ax.set_xticks(xlocs)
    _ = ax.set_xticklabels(xaxis, minor=True, fontsize=4)
    _ = ax.set_yticks(ylocs)
    _ = ax.set_yticklabels(yaxis, minor=True,fontsize=4)
    ax.set_title((fname + " conformations").upper())
    # save as pdf
    file_path = pathlib.Path(save_path) / f"{fname}_percentages.pdf"
    fig.savefig(fname=file_path, format="pdf")


def correlation_plot(data, fname, save_path, base="T", label_offset=0.05):
    # define colormap
    cmap = mpl.colors.ListedColormap([
        "blue",
        "cornflowerblue",
        "lightskyblue",
        "white",
        "mistyrose",
        "tomato",
        "red"])
    bounds = [-1.0, -.73, -.53, -.3, .3, .53, .73, 1.0]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cmap.set_bad(color="gainsboro")

    # reorder dataset
    coordinates = list(set(data.index.get_level_values(0)))
    data = data.loc[coordinates][coordinates].sort_index(
        level=1, axis=0).sort_index(level=1, axis=1)

    for crd1, crd2 in itertools.combinations_with_replacement(
            coordinates,
            r=2):
        crd_data = data.loc[crd1][crd2]

        # plot
        fig, ax = plt.subplots(
            1,
            1,
            dpi=300,
            tight_layout=True)
        im = ax.imshow(crd_data, cmap=cmap, norm=norm, aspect='auto')
        plt.colorbar(im)

        # axes
        units = set(crd_data.index)
        xlocs = np.arange(len(units))
        _ = ax.set_xticks(xlocs)
        _ = ax.set_xticklabels(units, rotation=90)

        ylocs = np.arange(len(units))
        _ = ax.set_yticks(ylocs)
        _ = ax.set_yticklabels(units)
        ax.set_title(f"rows: {crd1} | columns: {crd2}")
        plt.tight_layout()

        # save as pdf
        file_path = pathlib.Path(save_path) / f"{crd1}_{crd2}.pdf"
        fig.savefig(fname=file_path, format="pdf")

        plt.close()

    # plot
    fig, ax = plt.subplots(
        1,
        1,
        dpi=300,
        tight_layout=True)
    im = ax.imshow(data, cmap=cmap, norm=norm, aspect='auto')
    plt.colorbar(im)

    # axes
    start = len(data) // (2 * len(coordinates))
    step = 2 * start
    locs = np.arange(start, len(data)-1, step)
    _ = ax.set_xticks(locs)
    _ = ax.set_yticks(locs)
    _ = ax.set_xticklabels(coordinates, rotation=90)
    _ = ax.set_yticklabels(coordinates)

    plt.tight_layout()

    # save as pdf
    file_path = pathlib.Path(save_path) / f"{fname}.pdf"
    fig.savefig(fname=file_path, format="pdf")


def basepair_plot(
        data,
        fname,
        save_path,
        base="T",
        label_offset=0.05):
    # define colormap
    cmap = mpl.colors.ListedColormap([
        "blue",
        "cornflowerblue",
        "lightskyblue",
        "white",
        "mistyrose",
        "tomato",
        "red"])
    bounds = [-.6, -.5, -.4, -.3, .3, .4, .5, .6]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cmap.set_bad(color="gainsboro")

    category = data.index.to_series().apply(lambda s: s[0:4])
    data["category"] = category
    
    for cat in category.unique():
        cat_df = data[data["category"] == cat]
        cat_df = cat_df.drop("category", axis=1)
        # plot
        fig, ax = plt.subplots(
            1,
            1,
            dpi=300,
            figsize=(8, 7),
            tight_layout=True)
        im = ax.imshow(cat_df, cmap=cmap, norm=norm, aspect='auto')
        plt.colorbar(im)

        # axes
        xlocs = np.arange(len(cat_df.columns))
        _ = ax.set_xticks(xlocs)
        _ = ax.set_xticklabels(cat_df.columns.to_list(), rotation=90)

        ylocs = np.arange(len(cat_df.index))
        _ = ax.set_yticks(ylocs)
        _ = ax.set_yticklabels(cat_df.index.to_list(),fontsize=4)
        
        ax.set_title(
            f"Correlation for basepair group {cat}")
        plt.tight_layout()

        # save as pdf
        file_path = pathlib.Path(save_path) / f"{cat}.pdf"
        fig.savefig(fname=file_path, format="pdf")

        plt.close()

    data = data.sort_values(by="category")
    # cat_count = category.value_counts()
    # category = category.unique()
    # category.sort()
    data = data.drop("category", axis=1)

    # plot
    fig, ax = plt.subplots(
        1,
        1,
        dpi=300,
        figsize=(7.5, 5),
        tight_layout=True)
    im = ax.imshow(data, cmap=cmap, norm=norm, aspect='auto')
    plt.colorbar(im)

    # axes
    xlocs = np.arange(len(data.columns))
    _ = ax.set_xticks(xlocs)
    _ = ax.set_xticklabels(data.columns.to_list(), rotation=90)

    # if yaxis:
    # y_positions = [cat_count[category[i]] for i in range(len(category))]
    # ylocs = np.cumsum(y_positions)
    # _ = ax.set_yticks(ylocs)
    # _ = ax.set_yticklabels(category)
    # else:
    # _ = ax.set_yticklabels([])

    ax.set_title("Correlation for all basepairs")
    plt.tight_layout()

    # save as pdf
    file_path = pathlib.Path(save_path) / f"{fname}.pdf"
    fig.savefig(fname=file_path, format="pdf")
