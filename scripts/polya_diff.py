#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
import pandas as pd
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import seaborn as sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Overview report of poly(A) tail lengths.')
parser.add_argument(
    '-i', metavar='input', type=str, help="Input TSV.", required=True)
parser.add_argument(
    '-g', metavar='global', type=str, help="Global results TSV.", required=True)
parser.add_argument(
    '-t', metavar='per_tr', type=str, help="Per transcript results TSV.", required=True)
parser.add_argument(
    '-c', metavar='min_cov', type=int, help="Minimum coverage.", required=True)
parser.add_argument(
    '-r', metavar='report', type=str, help="Report PDF.", required=True)
parser.add_argument('-x', action="store_true",
                    help="Plot per-transcript distributions.", default=False)


def _make_distplot(data, title, label, xlab, ylab, pages):
    """ Make distplot with median. """
    ax = sns.distplot(data, kde=False, hist_kws={
                      "label": label}, norm_hist=False)
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_xlabel(xlab)
    ax.legend(loc='best')
    pages.savefig()
    plt.clf()


def _make_boxplot(df, by, title, pages, showfliers=False, mdf=None):
    """ Make box plot. """
    ax = sns.boxplot(x='sample', y='polya_length', hue=by,
                     showfliers=showfliers, data=df)
    ax.set_title(title)
    pages.savefig()
    plt.clf()


def _make_boxplot2(df, title, pages, showfliers=False, mdf=None):
    """ Make box plot2. """
    ax = sns.boxplot(x='group', y='polya_length',
                     showfliers=showfliers, data=df)
    mdf_d = OrderedDict()

    for r in mdf.itertuples():
        mdf_d[r.Index] = r.polya_length

    for group, med in mdf_d.items():
        ax.text(list(mdf_d.keys()).index(group), med + 0.5, np.round(med, 2),
                horizontalalignment='center', size='x-small', color='w', weight='semibold')
    ax.set_title(title)
    pages.savefig()
    plt.clf()


def _corrfunc(x, y, **kws):
    r, p = stats.pearsonr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.3f} p={:.3f}".format(r, p),
                xy=(.1, .9), xycoords=ax.transAxes)


def _make_pairplots(df, title, pages):
    """ Make pair plots. """
    controls = df[df.group == 'control']['sample'].unique()
    treatments = df[df.group == 'treatment']['sample'].unique()
    v = df[['contig', 'sample', 'polya_length']]
    v = v.reset_index().pivot_table(
        index=['contig'], columns='sample', values='polya_length').reset_index()
    v = v.groupby('contig').median()
    sns.pairplot(v[controls].dropna(), plot_kws={
                 's': 5, 'alpha': 0.65}).map_lower(_corrfunc)
    plt.subplots_adjust(top=0.9)
    plt.suptitle("Per-transcript tail length medians: control.")
    pages.savefig()
    plt.clf()
    sns.pairplot(v[treatments].dropna(), plot_kws={
                 's': 5, 'alpha': 0.65}).map_lower(_corrfunc)
    plt.suptitle("Per-transcript tail length medians: treatment.")
    plt.subplots_adjust(top=0.9)
    plt.suptitle(title)
    pages.savefig()
    plt.clf()


def _filter_medians(mdf, min_cov, pages):
    total_tr = len(mdf)
    pos = (mdf['count']['control'] >= min_cov) & (
        mdf['count']['treatment'] > min_cov)
    pass_tr = sum(pos)
    mdf_flt = mdf[pos]
    d = pd.DataFrame(OrderedDict(
        [('Category', ["Pass", "Fail"]), ('Count', [pass_tr, total_tr - pass_tr])]))
    d.plot.bar(x='Category', y='Count',
               title="Filtering for transcripts with counts >= {}".format(min_cov))
    plt.tight_layout()
    pages.savefig()
    plt.clf()
    return mdf_flt.copy()


if __name__ == '__main__':
    args = parser.parse_args()
    pages = PdfPages(args.r)
    sns.set_style("whitegrid")

    tails = pd.read_csv(args.i, sep="\t")
    mdf = tails[['contig', 'polya_length', 'group']].groupby(
        ['contig', 'group']).agg(['median', 'count'])
    mdf = mdf.unstack(level=-1)
    mdf.columns = mdf.columns.droplevel()

    mdf_flt = _filter_medians(mdf, args.c, pages)
    mdf_flt.loc[:, ('median', 'diff')] = mdf_flt.loc[
                    :, ('median', 'treatment')].values - mdf_flt.loc[:, ('median', 'control')].values
    mdf_flt = mdf_flt.sort_index(axis=1).copy()

    trs = mdf_flt.index.unique()

    tails_flt = tails[tails.contig.isin(trs)]
    _make_boxplot(tails[['contig', 'group', 'sample', 'polya_length']],
                  by='group', title="Boxplot overview (no outliers).", pages=pages)
    _make_boxplot(tails_flt[['contig', 'group', 'sample', 'polya_length']], by='group',
                  title="Boxplot overview after filtering (no outliers).", pages=pages)
    _make_pairplots(tails_flt[['contig', 'group', 'sample', 'polya_length']],
                    title="Data overview after filtering.", pages=pages)

    gu, gp = stats.mannwhitneyu(tails_flt[tails_flt.group == "control"].polya_length.values, tails_flt[
                                tails_flt.group == "treatment"].polya_length.values, alternative="two-sided")
    title = "Global dist. (no outliers): p-value={:.3f} U={:.3f}".format(
        gp, gu)
    tails_flt_med = tails_flt[
        ['group', 'polya_length']].groupby('group').median()

    _make_boxplot2(df=tails_flt, title=title, pages=pages,
                   showfliers=False, mdf=tails_flt_med)
    pd.DataFrame({'control': tails_flt_med.loc['control'].polya_length, 'treatment': tails_flt_med.loc['treatment'].polya_length, 'diff': tails_flt_med.loc[
                 'treatment'].polya_length - tails_flt_med.loc['control'].polya_length, 'p-value': [gp], 'U': [gu]}).to_csv(args.g, sep="\t", index=False)

    ulist, plist = [], []
    for tr in trs:
        trd = tails_flt[tails_flt.contig == tr]
        u, p = stats.mannwhitneyu(trd[trd.group == "control"].polya_length.values, trd[
                                  trd.group == "treatment"].polya_length.values, alternative="two-sided")
        ulist.append(u)
        plist.append(p)

    mdf_flt['u'] = ulist
    mdf_flt['p-value'] = plist
    mdf_flt['FDR'] = multipletests(plist, method="fdr_bh")[1]
    mdf_flt = mdf_flt.sort_values(by='FDR')
    mdf_flt.to_csv(args.t, sep="\t", index=True)

    _make_distplot(mdf_flt[('median', 'diff')].values, title="Distribution of median differences", label=mdf_flt[
                   ('median', 'diff')].median(), xlab="Tail length difference (treatment - control)", ylab="Count", pages=pages)
    _make_distplot(mdf_flt[mdf_flt.FDR < 0.05][('median', 'diff')].values, title="Distribution of median differences: FDR < 0.05", label=mdf_flt[
                   mdf_flt.FDR < 0.05][('median', 'diff')].median(), xlab="Tail length difference (treatment - control)", ylab="Count", pages=pages)

    if args.x:
        for tr in mdf_flt.index:
            tr_data = tails_flt[['contig', 'group', 'polya_length']].set_index('contig').loc[
                                                                               tr]
            title = "{} (no outliers): FDR={:.3f} U={:.3f}".format(
                tr, mdf_flt.loc[tr]['FDR'].values[0], mdf_flt.loc[tr]['u'].values[0])
            _make_boxplot2(df=tr_data, title=title, pages=pages,
                           showfliers=False, mdf=tr_data.groupby('group').median())

    pages.close()
