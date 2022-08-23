import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy.stats import gaussian_kde
import hddm
import os
from tqdm import tqdm

os.chdir('/Users/sean/documents/PICCOA/Models/HDDM')

plt.style.use('grayscale')
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['figure.figsize'] = (10.4, 6.8)
matplotlib.rcParams['lines.linewidth'] = 2.5
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=['0', '0.5']) 

pd.options.mode.chained_assignment = None  # ignore chaining warning; default='warn'               

# Load and clean ----------
dat = pd.read_csv('clean_dots.csv')
dat['subj_idx'] = dat['id']
dat['rt'] = dat['RT']  

dat['condition_c'] = np.where(dat.condition=='Stable', -.5, .5)
dat['age_group_c'] = np.where(dat.age_group=='YA', -.5, .5)
dat['trial0_c']    = dat.trial0 - .5
dat['colour0_c']   = dat.colour0 - .5

# Fit hddms ----------
regs = ['a ~ age_group_c * condition_c * trial0_c * colour0_c', 
        'v ~ age_group_c * condition_c * trial0_c * colour0_c', 
        't ~ age_group_c * condition_c * trial0_c * colour0_c']

ddm = hddm.models.HDDMRegressor(dat, regs)
ddm.find_starting_values()
ddm.sample(2000, burn=500, dbname=f'fits/traces_fullmod', db='pickle')
ddm.save(f'fits/fit_fullmod')

# load model, if already fit
ddm = hddm.load('fits/fit_fullmod')

# add p-values and gewekes
ddm_summary = ddm.gen_stats()

fig    = plt.figure(1, figsize=[30,30], tight_layout=True)
pars   = ddm_summary[~ddm_summary.index.str.contains('Intercept_')].index.to_list()
traces = ddm.nodes_db.node[[i for i in pars]]

ddm_summary['p_err']  = ddm_summary['mc err']/ddm_summary['std']
ddm_summary['geweke'] = np.nan
ddm_summary['P']      = np.nan

for i,p in enumerate(traces):
    trace = p.trace()

    # ax = fig.add_subplot(round(len(pars)/9), 9, i+1)
    # ax.plot(trace)
    # ax.set(xlabel='', ylabel='', title=pars[i])

    # geweke
    m1,m2 = np.mean(trace[:int(trace.shape[0]*.5)]), np.mean(trace[int(trace.shape[0]*.9):])
    v1,v2 = np.var(trace[:int(trace.shape[0]*.5)]), np.var(trace[int(trace.shape[0]*.9):])
    z = (m2-m1)/np.sqrt(v2+v1)

    # p-value
    p_lt_0 = np.mean(trace < 0)
    p      = p_lt_0 if ddm_summary.loc[pars[i]]['mean'] < 0 else 1-p_lt_0

    # save
    ddm_summary.loc[pars[i], 'geweke'] = z
    ddm_summary.loc[pars[i], 'P'] = 1-p # to match frequentist interpretation

# fig.savefig('figs/traceplots_fullfit.png')
ddm_summary.to_csv('stats/summary_fullmod_all.csv')
ddm_summary[~ddm_summary.index.str.contains('Intercept_')].to_csv('stats/summary_fullmod_main.csv')

# Error bar plot
titles = {'a':'Decision Threshold', 'v':'Drift Rate', 't':'Non-Decision Time'}

fig,ax = plt.subplots(1,3, figsize=[25,10])

for i,p in enumerate(['a','v','t']):
    tmp = ddm_summary[ddm_summary.index.str[0]==p]
    tmp = tmp[~tmp.index.str.contains('Intercept')]
       
    tmp['lb'] = tmp['mean'] - tmp['2.5q']
    tmp['ub'] = tmp['97.5q'] - tmp['mean']
    ci = tmp[['lb','ub']].T.to_numpy()

    ax[i].errorbar(tmp.index.to_list(), tmp['mean'], fmt='o', yerr=ci, ms=10)
    ax[i].axhline(0, linestyle='--', color='grey', alpha=0.7)
    ax[i].set_xticklabels(tmp.index.to_list(), rotation=90, size=12)
    ax[i].tick_params(axis='y', which='major', labelsize=12)
    ax[i].set(title=titles[p])

    if i ==0:
        ax[i].set(ylabel=r'$\beta$')

fig.subplots_adjust(bottom=.4)

plt.show()
fig.savefig('figs/summary_point_all.png')
plt.close()

# Posterior predictive check ---------
X = {'age_group':[], 'condition':[], 'trial':[], 'colour':[]}

for a in ['YA','OA']:
    for c in ['Stable','Decreasing']:
        for t in np.linspace(0, 1, num=800):
            for col in np.linspace(0, 1, num=100):
                X['age_group'].append(a)
                X['condition'].append(c)
                X['trial'].append(t)
                X['colour'].append(col)

X = pd.DataFrame(X)
X['condition_c'] = np.where(X.condition=='Stable', -.5, .5)
X['age_group_c'] = np.where(X.age_group=='YA', -.5, .5)
X['trial0_c']    = X.trial - .5
X['colour0_c']   = X.colour - .5

for p in ['a','v','t']:
    tmp = ddm_summary[ddm_summary.index.str[0]==p]
    tmp = tmp[~tmp.index.str.contains('Intercept')]

    pred = np.zeros((X.shape[0], tmp.shape[0]))
    for i,par in enumerate(tmp.index.to_list()):
        b = tmp.loc[par,'mean']
        pname = par[2:]

        if ':' in pname:
            varnames = pname.split(':')
            out = X[varnames[0]].copy()
            for v in varnames[1:]:
                out*=X[v]

            X[pname] = out

        pred[:,i] = b*X[pname]

    b0   = ddm_summary.loc[f'{p}_Intercept','mean']
    X[p] = b0 + pred.sum(axis=1)
    X[f'{p}_z'] = (X[p]-X[p].mean())/X[p].std()


# Parameter predictions (age_group)--------
titles = {'a':'Decision Threshold', 'v':'Drift Rate', 't':'Non-Decision Time'}
fig,ax = plt.subplots(1,3,figsize=[15,5])

for i,p in enumerate(['a','v','t']):

    y = X[p]*1000 if p=='t' else X[p]
    yrange = y.min()-.2, y.max()+.2

    sns.pointplot(x=X['age_group'], y=y, hue=X['condition'], 
                  linestyles='', ci='sd', ax=ax[i], 
                  palette=('darkred','orange'))

    ax[i].set(xlabel='', ylabel='', title=titles[p], ylim=yrange, 
              xticks=[0,1], xticklabels=['Young','Old'])

    if i==0:
        ax[i].legend(title='Condition')
    else:
        ax[i].legend([],[], frameon=False)

plt.show()
fig.savefig('figs/parameter_predictions.png')
plt.close()


# Continuous plots
X['stimbin'] = pd.qcut(X['colour'], 10, labels=False)
X['timebin'] = pd.qcut(X['trial'], 4, labels=False)

titles  = {'a':'Decision Threshold', 'v':'Drift Rate', 't':'Non-Decision Time'}
yranges = {'a': (X.a.min(), X.a.max()), 'v':(X.v.min(), X.v.max()), 't':(X.t.min(), X.t.max())}
ag_lab  = {'YA':'Young', 'OA':'Old'}

for par in ['a','v','t']:
    fig     = plt.figure(1, figsize=[10,6], tight_layout=True)
    counter = 1
    for ag in ('YA', 'OA'):   
        for cond in ('Stable', 'Decreasing'):
            tmp = X[(X.timebin.isin([0,3])) & (X.condition==cond) * (X.age_group==ag)]
            ax  = fig.add_subplot(2, 2, counter)

            sns.lineplot(x=tmp['stimbin'], y=tmp[par], hue=tmp['timebin'], palette=['red','blue'], ci='sd', ax=ax)
            ax.set(xlabel='', ylabel='', title=f'{ag_lab[ag]}\n{cond}', 
                xticks=[0,9], xticklabels=['Very Purple', 'Very Blue'], 
                ylim=yranges[par])

            if par=='v':
                ax.axhline(0, linestyle='--', color='grey', alpha=0.5)

            if counter==1:
                ax.legend(['First 200 trials', 'Last 200 trials'], prop={'size': 10})
            else:
                ax.legend([],[], frameon=False)

            counter+=1

    fig.suptitle(titles[par])
    plt.show()
    fig.savefig(f'figs/cont_plot_{par}.png')
    plt.close()


