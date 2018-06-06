import os
import numpy as np
import pandas as pd
import peakutils
from glob import glob


def get_df(filetype, filename, special=None):
    if filetype == 'emc2':
        if special != 'cm':
            names = ['T', 'mu', 'E-mu*x', 'x', 'phi', 'var(E)', 'var(x)',
                     'E_lte-mu*x_lte', 'x_lte', 'phi_lte',
                     'E_mf-mu*x_mf', 'x_mf', 'phi_mf',
                     'E_hte-mu*x_hte', 'x_hte', 'phi_hte', 'lro']
            try:
                df = pd.read_table(filename, names=names, usecols=range(17), sep='\s+', header=None, index_col=['T', 'mu'], dtype=float)
            except AttributeError:
                df = pd.DataFrame(columns=names).set_index(['T', 'mu'])

            # df['E'] = df['E-mu*x'] + df.index.levels[1][df.index.labels[1]] * df['x']
            # df['G'] = df['phi'] + df.index.levels[1][df.index.labels[1]] * df['x']
            if special == 'innerT':
                df = df.swaplevel(0, 1)
        else:
            names = ['T', 'mu', 'E', 'x', 'G', 'var(E)', 'var(x)',
                     'E_lte', 'x_lte', 'G_lte',
                     'E_mf', 'x_mf', 'G_mf',
                     'E_hte', 'x_hte', 'G_hte', 'lro']
            names_seletecd = ['E', 'G', 'var(E)', 'lro']
            try:
                df = pd.read_table(filename, names=names, usecols=range(17), sep='\s+', header=None, index_col=['x', 'T'], dtype=float)[names_seletecd]
            except AttributeError:
                df = pd.DataFrame(columns=names).set_index(['x', 'T'])[names_seletecd]

    elif filetype == 'phb':
        names = ['T', 'mu', 'x1', 'x2', 'E1', 'E2']
        try:
            df = pd.read_table(filename, names=names, sep='\s+', header=None, index_col=['T'])
        except KeyError:
            df = pd.read_table(filename, names=names + ['lro1', 'lro2'], sep='\s+', header=None, index_col=['T'])

    elif filetype in ['gs', 'fit', 'predstr']:
        if filetype == 'gs':
            if special != 'mmaps':
                names = ['c', 'E', 'E_fit', 'str_id']
                df = pd.read_table(filename, names=names, sep='\s+', header=None)
            else:
                df = pd.read_table(filename, sep='\s+', header=None)
                num_cols = len(df.columns)
                names = ['c' + str(i) for i in range(num_cols - 4)] + ['E', 'E_fit', 'E_diff', 'str_id']
                df.columns = names

        elif filetype == 'fit':
            if special != 'mmaps':
                names = ['c', 'E', 'E_fit', 'E_diff', 'wt', 'str_id']
                df = pd.read_table(filename, names=names, sep='\s+', header=None)
            else:
                df = pd.read_table(filename, sep='\s+', header=None)
                num_cols = len(df.columns)
                names = ['c' + str(i) for i in range(num_cols - 5)] + ['E', 'E_fit', 'E_diff', 'wt', 'str_id']
                df.columns = names

        elif filetype == 'predstr':
            if special != 'mmaps':
                names = ['c', 'E', 'E_fit', 'str_id', 'status']
                df = pd.read_table(filename, names=names, sep='\s+', header=None)
            else:
                df = pd.read_table(filename, sep='\s+', header=None)
                num_cols = len(df.columns)
                names = ['c' + str(i) for i in range(num_cols - 3)] + ['E_fit', 'str_id', 'status']
                df.columns = names

        df.index.name = '#'
        df = df.reset_index().set_index('str_id')

        if 'E' in df.columns:
            if 'E_diff' not in df.columns:
                df['E_diff'] = df['E'] - df['E_fit']
            df['E_diff_abs'] = df['E_diff'].abs()

        if filetype == 'gs':
            dirname = os.path.dirname(filename)
            filename_ref = os.path.join(dirname, 'ref_energy.out')
            if os.path.isfile(filename_ref):
                ref = np.loadtxt(filename_ref)
                if special != 'mmaps':
                    addition = (1-df['c']) * ref[0] + df['c'] * ref[1]
                else:
                    addition = 0
                    for i in [df['c' + str(i)] * ref[i] for i in range(len(ref))]:
                        addition += i
                if 'E' in df.columns:
                    df['E_ref_0'] = df['E'] + addition
                df['E_fit_ref_0'] = df['E_fit'] + addition

    elif filetype == 'clusinfo':
        names = ['n_pts', 'd', 'multi', 'eci']
        df = pd.read_table(filename, names=names, sep='\s+', header=None)

    return df


def get_peak_idx_array(var_diff_abs, thres, min_dist, thres_abs):
    # var_diff_abs -= var_diff_abs[var_diff_abs < var_diff_abs.quantile(0.8)].mean()
    var_diff_abs_mean = var_diff_abs[var_diff_abs < var_diff_abs.quantile(0.8)].mean()
    argmax = np.argmax(var_diff_abs.reset_index(drop=True))
    argmax = argmax if not np.isnan(argmax) and var_diff_abs.iloc[argmax] > thres_abs else None
    peak_idx_array = peakutils.indexes(var_diff_abs, thres, min_dist)
    mask = (((var_diff_abs.iloc[peak_idx_array] - var_diff_abs_mean)/var_diff_abs.max() > thres) \
           & (var_diff_abs.iloc[peak_idx_array] > thres_abs)).values
    peak_idx_array = peak_idx_array[mask]
    return peak_idx_array, argmax


def get_boundary_idx(var_diff_abs, peak_idx_pair, bd_tol):
    var_diff_abs = var_diff_abs - var_diff_abs[var_diff_abs < var_diff_abs.quantile(0.8)].mean()
    boundary_idx_pair = [peak_idx_pair[0] - 1, peak_idx_pair[1]]
    i = -1
    while peak_idx_pair[0] + i > 0 and var_diff_abs.iloc[peak_idx_pair[0] + i] > bd_tol * var_diff_abs.iloc[peak_idx_pair[0]]:
        boundary_idx_pair[0] = peak_idx_pair[0] + i - 1
        i -= 1
    i = 1
    while peak_idx_pair[1] + i < len(var_diff_abs) and var_diff_abs.iloc[peak_idx_pair[1] + i] > bd_tol * var_diff_abs.iloc[peak_idx_pair[1]]:
        boundary_idx_pair[1] = peak_idx_pair[1] + i
        i += 1
    return boundary_idx_pair


def plot_emc2(df, trans_df, columns):
    cycle = df.index.levels[0][np.unique(df.index.labels[0])]
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, len(columns), sharex=True, squeeze=False, figsize=[3.5 * len(columns), 6])
    xaxis = df.index.names[1]
    for idx, col in enumerate(columns):
        axes[0, idx].set_title(col)
        axes[1, idx].set_title(col + '_diff')
        axes[1, idx].set_xlabel(xaxis)

    # colors = plt.rcParams['axes.color_cycle']
    colors = [i['color'] for i in plt.rcParams['axes.prop_cycle'].__dict__['_left']]
    for idx_cyc, cyc in enumerate(cycle):
        df_at_cyc = df.xs(cyc)
        trans_df_at_cyc = trans_df.xs(cyc) if cyc in trans_df.index else None
        trans_df_at_cyc = trans_df_at_cyc.to_frame().T if not isinstance(trans_df_at_cyc, pd.DataFrame) and not trans_df_at_cyc is None else trans_df_at_cyc
        color = colors[idx_cyc % len(colors)]
        for idx, col in enumerate(columns):
            plt.sca(axes[0, idx])
            plt.plot(df_at_cyc.index, df_at_cyc[col], '-x', ms=3, color=color, label=df.index.names[0] + ' = ' + str(cyc))
            if not trans_df_at_cyc is None:
                plt.plot(trans_df_at_cyc[xaxis + '1'], trans_df_at_cyc[col + '1'].values, 'o', color=color, ms=4)
                plt.plot(trans_df_at_cyc[xaxis + '2'], trans_df_at_cyc[col + '2'].values, 's', color=color, ms=4)

            plt.sca(axes[1, idx])
            plt.plot(df_at_cyc.index, df_at_cyc[col].diff(), '-x', ms=3, color=color)
            if not trans_df_at_cyc is None:
                plt.plot(df_at_cyc.index[trans_df_at_cyc['pk_idx1'].astype(int)], df_at_cyc[col].diff().iloc[trans_df_at_cyc['pk_idx1'].astype(int)].values, 'o', color=color, ms=4)
                plt.plot(df_at_cyc.index[trans_df_at_cyc['pk_idx2'].astype(int)], df_at_cyc[col].diff().iloc[trans_df_at_cyc['pk_idx2'].astype(int)].values, 's', color=color, ms=4)

    legend_ncol = 2 if len(cycle) > 10 else 1
    axes[0, -1].legend(loc=0, fontsize='x-small', handlelength=1, labelspacing=0.1, ncol=legend_ncol)
    plt.autoscale(tight=True)
    plt.locator_params(nbins=5)
    plt.tight_layout(h_pad=0.1, w_pad=0.1)


def get_emc2_transition(df, var='lro', thres=0.6, min_dist=2, thres_abs=0.1, bd_tol=0.3, return_single=False, cm=False, plot=False):
    cycle = df.index.levels[0][np.unique(df.index.labels[0])]
    outer = df.index.names[0]
    inner = df.index.names[1]
    trans_list = []
    for cyc in cycle:
        df_at_cyc = df.xs(cyc)
        var_diff_abs = df_at_cyc[var].diff().abs()
        peak_idx_array, argmax = get_peak_idx_array(var_diff_abs, thres, min_dist, thres_abs)
        # get the first and last peak points to characterize the transition range
        if len(peak_idx_array):
            if not return_single:
                for peak_idx in peak_idx_array:
                    peak_idx_pair = [peak_idx, peak_idx]
                    boundary_idx_pair = get_boundary_idx(var_diff_abs, peak_idx_pair, bd_tol)
                    trans_dict = {outer: cyc,
                                  inner + '1': df_at_cyc.index[boundary_idx_pair[0]],
                                  inner + '2': df_at_cyc.index[boundary_idx_pair[1]],
                                  'lro1': df_at_cyc['lro'].iloc[boundary_idx_pair[0]],
                                  'lro2': df_at_cyc['lro'].iloc[boundary_idx_pair[1]],
                                  'bd_idx1': boundary_idx_pair[0], 'bd_idx2': boundary_idx_pair[1],
                                  'pk_idx1': peak_idx_pair[0], 'pk_idx2': peak_idx_pair[1],
                                  'argmax': argmax}
                    if not cm:
                        trans_dict.update({'x1': df_at_cyc['x'].iloc[boundary_idx_pair[0]], 'x2': df_at_cyc['x'].iloc[boundary_idx_pair[1]],
                            'E-mu*x1': df_at_cyc['E-mu*x'].iloc[boundary_idx_pair[0]], 'E-mu*x2': df_at_cyc['E-mu*x'].iloc[boundary_idx_pair[1]],
                            'phi1': df_at_cyc['phi'].iloc[boundary_idx_pair[0]], 'phi2': df_at_cyc['phi'].iloc[boundary_idx_pair[1]]})
                    else:
                        trans_dict.update({'E1': df_at_cyc['E'].iloc[boundary_idx_pair[0]], 'E2': df_at_cyc['E'].iloc[boundary_idx_pair[1]],
                            'G1': df_at_cyc['G'].iloc[boundary_idx_pair[0]], 'G2': df_at_cyc['G'].iloc[boundary_idx_pair[1]]})
                    trans_list.append(trans_dict)
            else:
                # find the peak width
                peak_idx_pair = peak_idx_array[[0, -1]]
                boundary_idx_pair = get_boundary_idx(var_diff_abs, peak_idx_pair, bd_tol)
                trans_dict = {outer: cyc,
                              inner + '1': df_at_cyc.index[boundary_idx_pair[0]],
                              inner + '2': df_at_cyc.index[boundary_idx_pair[1]],
                              'lro1': df_at_cyc['lro'].iloc[boundary_idx_pair[0]],
                              'lro2': df_at_cyc['lro'].iloc[boundary_idx_pair[1]],
                              'bd_idx1': boundary_idx_pair[0], 'bd_idx2': boundary_idx_pair[1],
                              'pk_idx1': peak_idx_pair[0], 'pk_idx2': peak_idx_pair[1],
                              'argmax': argmax}
                if not cm:
                    trans_dict.update({'x1': df_at_cyc['x'].iloc[boundary_idx_pair[0]], 'x2': df_at_cyc['x'].iloc[boundary_idx_pair[1]],
                        'E-mu*x1': df_at_cyc['E-mu*x'].iloc[boundary_idx_pair[0]], 'E-mu*x2': df_at_cyc['E-mu*x'].iloc[boundary_idx_pair[1]],
                        'phi1': df_at_cyc['phi'].iloc[boundary_idx_pair[0]], 'phi2': df_at_cyc['phi'].iloc[boundary_idx_pair[1]]})
                else:
                    trans_dict.update({'E1': df_at_cyc['E'].iloc[boundary_idx_pair[0]], 'E2': df_at_cyc['E'].iloc[boundary_idx_pair[1]],
                        'G1': df_at_cyc['G'].iloc[boundary_idx_pair[0]], 'G2': df_at_cyc['G'].iloc[boundary_idx_pair[1]]})
                trans_list.append(trans_dict)


    if not cm:
        trans_df = pd.DataFrame(trans_list, columns=[outer, inner + '1', inner + '2', 'x1', 'x2', 'lro1', 'lro2', 'E-mu*x1', 'E-mu*x2', 'phi1', 'phi2',
                         'bd_idx1', 'bd_idx2', 'pk_idx1', 'pk_idx2', 'argmax']).set_index(outer)
    else:
        trans_df = pd.DataFrame(trans_list, columns=[outer, inner + '1', inner + '2', 'lro1', 'lro2', 'E1', 'E2', 'G1', 'G2',
                             'bd_idx1', 'bd_idx2', 'pk_idx1', 'pk_idx2', 'argmax']).set_index(outer)

    if plot:
        if not cm:
            plot_emc2(df, trans_df, columns=['lro', 'x'])
        else:
            plot_emc2(df, trans_df, columns=['lro'])

    return trans_df


def get_emc2_last_points(df, cm=False):
    cycle = df.index.levels[0][np.unique(df.index.labels[0])]
    outer = df.index.names[0]
    inner = df.index.names[1]
    trans_list = []
    for cyc in cycle:
        df_at_cyc = df.xs(cyc)
        trans_dict = {outer: cyc,
                      inner + '1': df_at_cyc.index[-2],
                      inner + '2': df_at_cyc.index[-1],
                      'lro1': df_at_cyc['lro'].iloc[-2],
                      'lro2': df_at_cyc['lro'].iloc[-1]}
        if not cm:
            trans_dict.update({'x1': df_at_cyc['x'].iloc[-2], 'x2': df_at_cyc['x'].iloc[-1],
                'E-mu*x1': df_at_cyc['E-mu*x'].iloc[-2], 'E-mu*x2': df_at_cyc['E-mu*x'].iloc[-1],
                'phi1': df_at_cyc['phi'].iloc[-2], 'phi2': df_at_cyc['phi'].iloc[-1]})
        else:
            trans_dict.update({'E1': df_at_cyc['E'].iloc[-2], 'E2': df_at_cyc['E'].iloc[-1],
                'G1': df_at_cyc['G'].iloc[-2], 'G2': df_at_cyc['G'].iloc[-1]})
        trans_list.append(trans_dict)

    if not cm:
        trans_df = pd.DataFrame(trans_list, columns=[outer, inner + '1', inner + '2', 'x1', 'x2', 'lro1', 'lro2', 'E-mu*x1', 'E-mu*x2', 'phi1', 'phi2']).set_index(outer)
    else:
        trans_df = pd.DataFrame(trans_list, columns=[outer, inner + '1', inner + '2', 'lro1', 'lro2', 'E1', 'E2', 'G1', 'G2']).set_index(outer)

    return trans_df


def get_emc2_glob_df(glob_str, special=None):
    glob_df = pd.concat([get_df('emc2', filename, special=special) for filename in sorted(glob(glob_str))])
    return glob_df


def plot_phb_by_emc2(trans_df, special=None, kwargs={}, connect=False):
    import matplotlib.pyplot as plt
    cm = True if special == 'cm' else False
    if len(trans_df):
        kwargs1 = kwargs.copy()
        kwargs2 = kwargs1.copy()
        # if 'l' in glob_str or 'b' in glob_str:
            # kwargs1.update({'label': glob_str})
            # if connect:
                # kwargs1.update({'ls': '-'})
                # kwargs2.update({'dashes': [2, 2]})
        # elif 'r' in glob_str or 't' in glob_str:
            # kwargs2.update({'label': glob_str})
            # if connect:
                # kwargs1.update({'dashes': [2, 2]})
                # kwargs2.update({'ls': '-'})
        if connect:
            kwargs1.update({'ls': '-'})
            kwargs2.update({'dashes': [2, 2]})

        trans_df.reset_index(inplace=True)
        if not cm:
            if special != 'innerT':
                plt.plot(trans_df['x1'], trans_df['T'].values, 'x', ms=3, **kwargs1)
                plt.plot(trans_df['x2'], trans_df['T'].values, 'x', ms=3, **kwargs2)
            else:
                plt.plot(trans_df['x1'], trans_df['T1'].values, 'x', ms=3, **kwargs1)
                plt.plot(trans_df['x2'], trans_df['T2'].values, 'x', ms=3, **kwargs2)
        else:
            plt.plot(trans_df['x'], trans_df['T1'].values, 'x', ms=3, **kwargs1)
            plt.plot(trans_df['x'], trans_df['T2'].values, 'x', ms=3, **kwargs2)


def plot_phb(filelist, kwargs={}):
    import matplotlib.pyplot as plt
    colors = plt.rcParams['axes.color_cycle']
    kwargs_input = kwargs
    for idx, filename in enumerate(filelist):
        kwargs = {'color': colors[idx % len(colors)]} if not kwargs_input else kwargs_input
        phb_df = get_df('phb', filename)
        plt.plot(phb_df['x1'], phb_df.index.values, '-x', ms=3, label=filename, **kwargs)
        plt.plot(phb_df['x2'], phb_df.index.values, '-x', ms=3, **kwargs)
