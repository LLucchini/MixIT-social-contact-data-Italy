import os

import pandas as pd
import numpy as np
import itertools
import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.dates as mdates

import seaborn as sns
import scipy.stats as scs


font_dirs = ['/Library/Fonts/']
font_files = mpl.font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    mpl.font_manager.fontManager.addfont(font_file)
mpl.rcParams['font.family'] = 'Barlow'
sns.set(font='Barlow')


palette = color=['#223E5C','#5EACA3','#F1D670','#DC6046']
palette_cb = ['#364B9A','#6EA6CD','#C2E4EF','#EAECCC','#FEDA8B','#F67E4B','#A50026']
palette_cb_long = ['#364B9A','#4A7BB7','#6EA6CD','#98CAE1','#C2E4EF','#EAECCC','#FEDA8B','#FDB366','#F67E4B','#DD3D2D','#A50026']

location_map = {
    'conviventi':'home',
    'home':'home',
    'homeguest':'leisure',
    'work':'work',
    'school':'school',
    'leisure':'leisure',
    'shopping':'other',
    'restaurant':'leisure',
    'transport':'transport',
    'otherindoor':'other',
    'otheroutdoor':'other'
}

def set_plot_style(ax, xlim=None, ylim=None, panel_name=None, panel_loc=(-0.1,1.1), title=None, xlabel=None, ylabel=None, xcolor='#353535', ycolor='#353535', legend_frame = False, index_dates = False, legend=False, dates_formatter = 'months', dates_to_plot=None, freq=2,xlims=(None,None),offset_vbar = 0, date_names=False, date_names_offset=1.2, letters = [chr(i) for i in range(93,123)],fontsize=11,grid=False):
    if title:        ax.set_title(title,fontsize=fontsize+4)
    if not panel_name is None: ax.text(s=panel_name, x=panel_loc[0], y=panel_loc[1], rotation=0,  transform=ax.transAxes)
    if xlabel:       ax.set_xlabel(xlabel, fontsize=fontsize+2,fontweight='semibold',color=xcolor)
    if ylabel:       ax.set_ylabel(ylabel, fontsize=fontsize+2,fontweight='semibold',color=ycolor)
    if xlims:         ax.set_xlim(xlims)
    if xlim:         ax.set_xlim(xlim)
    if ylim:         ax.set_ylim(ylim)
    if index_dates:
        if dates_formatter == 'months':
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=freq))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        elif dates_formatter == 'days':
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=freq))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        elif dates_formatter == 'years':
            ax.xaxis.set_major_locator(mdates.YearLocator(interval=freq))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right',c=xcolor,fontsize=fontsize)
    ax.yaxis.set_tick_params(labelcolor=ycolor,labelsize=fontsize)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.minorticks_off()
    # add vertical lines for important dates
    if dates_to_plot:
        for i,(day, name) in enumerate(dates_to_plot.items()):
            if (pd.to_datetime(day,format='%d-%m-%Y',utc=True)>xlims[0])&(pd.to_datetime(day,format='%d-%m-%Y',utc=True)<=xlims[1]):
                ax.axvline(pd.to_datetime(day,format='%d-%m-%Y'), linestyle='--', color ='#757575', linewidth=1)
                if date_names: ax.text(s=name, x=pd.to_datetime(day,format='%d-%m-%Y')-np.timedelta64(offset_vbar,'D'), y=0.05, rotation=90, transform=ax.get_xaxis_transform(), fontsize=fontsize)
                else:
                    if letters != None:
                        ax.text(s=letters[i], x=pd.to_datetime(day,format='%d-%m-%Y')-np.timedelta64(offset_vbar,'D'), y=1.01, rotation=0, transform=ax.get_xaxis_transform(), fontsize=fontsize)
                        ax.text(s=letters[i]+': '+name, x=date_names_offset, y=1.5-0.1*(i+1), rotation=0, transform=ax.transAxes, fontsize=fontsize)
                    else: ax.text(s=[chr(i) for i in range(93,123)][i], x=pd.to_datetime(day,format='%d-%m-%Y')-np.timedelta64(offset_vbar,'D'), y=1.01, rotation=0, transform=ax.get_xaxis_transform(), fontsize=fontsize+2)
    if legend: ax.legend(frameon=legend_frame)
    else: ax.legend([],frameon=legend_frame)
    # remove grid
    ax.grid(grid)


##################################
def plot_istat_distributions(survey_data, istat_data, istat_col, ax=None, color_palette=palette_cb, subset_data=None, box_color='red', surve_istat_denoms = None, base_col='part_age_group_pop_5', outer_with_pop=False):
    s_data = survey_data.copy()
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,5))
    if type(istat_col) == str:
        istat_col = [istat_col]

    if not subset_data is None:
        s_data = subset_data.copy()
    # compute denominator for population fraction
    # survey denominator
    s_denom = (s_data[istat_col+['part_id']]
         .drop_duplicates(istat_col+['part_id'])
         .dropna()['part_id']
         .nunique()
        )
    # istat data denominator
    p = istat_data.groupby(istat_col, observed=False)[base_col].sum().to_frame()
    i_denom = p.sum()
    if not surve_istat_denoms is None:
        s_denom, i_denom = surve_istat_denoms
    r = ((s_data[istat_col+['part_id']]
         .drop_duplicates(istat_col+['part_id'])
         .dropna()
         .groupby(istat_col, observed=False)['part_id']
         .nunique()/ s_denom
        ).reset_index()
      )
    r.columns = istat_col+['survey data']

    respondent_dict = {f'{i}_x':f'respondent\n{iscol[i]}' for i, iscol in enumerate(istat_col)}
  
    p = p/i_denom

    if outer_with_pop: toBar = r.merge(p.dropna(), on=istat_col, how='outer').set_index(istat_col)[['survey data',base_col]]
    else:              toBar = r.merge(p.dropna(), on=istat_col).set_index(istat_col)[['survey data',base_col]]

    toBar[f'survey data'].plot(kind='bar', ax=ax, width=0.65, color=color_palette)#['#20175FFF','#82D8B0'])
    # toBar[f'survey data'].plot(kind='bar', ax=ax, width=0.65)#['#20175FFF','#82D8B0'])

    # ax_right = ax.twinx()
    ax.bar(range(len(toBar)), height=toBar[base_col], width=0.8, color='#00000000', edgecolor=box_color, label='national data')
    # ax.set_ylabel('average contacts', color='#20175FFF')
    ax.set_ylabel('fraction of population')
    ax.legend(frameon=False)
    # ax_right.legend(frameon=False, loc=(0.01,0.8))
    return ax

def zip_col_vals(x):
    nan_check = None
    for e in x.values:
        if pd.isna(e):      return np.nan
    if nan_check is None: return tuple(x.tolist())

#####################################################################
#####################################################################
######################## Symmetrization #############################
#####################################################################
def symmetrize_matrix(contacts,respondent_info,group_col,group_pop, order_col, sens=12, subset_group=None, subset_value=None, verbose=False, index_NA=None):
    if not subset_group is None:
        sub_respondents = respondent_info[respondent_info[subset_group]==subset_value]
        sub_contacts = contacts[(contacts['part_id'].isin(sub_respondents['part_id']))]
    else:
        sub_respondents = respondent_info.copy()
        sub_contacts = contacts[contacts['part_id'].isin(sub_respondents['part_id'])].copy()
    # get group index, sort it and use it to loc rows
    if index_NA is None:
        if group_col!=order_col:
            group_loc = contacts[[group_col,order_col,group_pop]].dropna().drop_duplicates().sort_values(order_col).loc[:,group_col].dropna().unique()
        else:
            group_loc = contacts[[group_col,group_pop]].dropna().drop_duplicates().sort_values(order_col).loc[:,group_col].dropna().unique()
        index_NA = contacts[[group_col,order_col,group_pop]].dropna().drop_duplicates().sort_values(order_col).loc[:,[group_col, group_pop]].dropna().drop_duplicates().set_index(group_col)[group_pop].to_dict()
    else:
        group_loc = list(index_NA.keys())
    # building total contact matrix
    Mij = sub_contacts.groupby(['part_id',group_col,f'cnt{group_col[4:]}'], as_index=False, dropna=False)['cont_id'].count()
    Mij = Mij.merge(pd.DataFrame(data=itertools.product(contacts['part_id'].unique(),contacts[f'cnt{group_col[4:]}'].dropna().unique()), 
                                 columns=['part_id',f'cnt{group_col[4:]}']), 
                    on=['part_id',f'cnt{group_col[4:]}'], how='right')
    Mij[group_col] = Mij['part_id'].map(contacts[['part_id',group_col]].dropna().drop_duplicates().set_index('part_id').to_dict()[group_col])
    Mij = Mij.dropna()
    Mij = Mij.groupby([group_col,f'cnt{group_col[4:]}'], as_index=False, dropna=False)['cont_id'].sum().pivot(index=group_col,columns=f'cnt{group_col[4:]}',values='cont_id')
    # fill missing indices or columns (reoports nans in output if) (T_{a,\tilde{a}})
    Mij = Mij.reindex(index=group_loc, columns=group_loc)

    # population group vector
    Ni = pd.Series(index_NA)
    # respondent_population group vector
    ni = sub_respondents.groupby([group_col], as_index=False)[group_pop].count().rename(columns={group_col:'group_col'}).set_index('group_col')[group_pop]
    ni = ni.reindex(index=group_loc).fillna(0)

    # Total contacts per group (C_{a,\tilde{a}})
    Mij = Mij.divide(ni, axis=0).fillna(0)
    # Symmetrized per capita contacts per group: matrix and df
    Mij_sym = ((Mij.T * ni * Ni).T + (Mij.T * ni * Ni)).replace(0,np.nan)
    PaPa_sym = pd.DataFrame((np.ones((Mij.shape)) * ni.values).transpose() + (np.ones((Mij.shape)) * ni.values), columns=ni.index, index=ni.index).replace(0,np.nan)
    # (T_{a,\tilde{a}}^{corrected})
    sym = (Mij_sym.divide(PaPa_sym))
    # get average contacts per individual after symmetrization
    df_sym = sym.divide(Ni, axis=0)
    df_sym.columns.name = 'contacts'
    df_sym.index.name = 'respondent'

    # Symmetrization check
    im, jm = sym.shape
    # sens = 12 # decimal unit to check for symmetry
    if verbose:
        sym_test = sym.values
        out = None
        for i in range(im):
            for j in range(jm):
                if round(sym_test[i,j],sens) != round(sym_test[j,i],sens): out = f'Matrix is not symmetric up to {sens} decimal number'
        if out is None: print(group_pop, subset_group, subset_value, 'Matrix is Symmetric!')
        else: print(out)
        
    np.seterr(over='ignore', invalid='ignore')
    Q = (np.diag(sym/sym.sum(0)).sum()-1)/(sym.shape[0]-1)
    
    return sym, df_sym.T.loc[group_loc[::-1],group_loc], Q # group_loc, Mij

#####################################################################
#####################             ##################################
def bootstrap_contact_matrices(rsi_all,respondentToPlot, contactToPlot, group_col, group_pop, order_col, group_loc=None, index_NA=None, bootstrap=100, bootstrap_dir='',save_bootstraps=True, subset_group=None, subset_value=None, sens=12, round_notes=1, minimal_output=True):
    if index_NA is None: index_NA = contactToPlot.drop_duplicates(group_col).sort_values(order_col).set_index(group_col)[group_pop].to_dict()
    if group_pop not in respondentToPlot.columns:
        respondentToPlot[group_pop] = respondentToPlot[group_col].map(contactToPlot[[group_col,group_pop]].drop_duplicates().set_index(group_col).to_dict()[group_pop])
    q1gs = []
    tmp1s = []
    Qs = []
    weights = []
    for i in range(bootstrap):
        print(f'Bootstrapping ... {i+1} / {bootstrap}\t', end='\r')
        rsi = respondentToPlot.sample(respondentToPlot.shape[0], replace=True).copy()
        rsi['new_id'] = rsi.reset_index(drop=True).index
        # modifying contacts to wellbehave under duplicare respondent sampling
        cts = contactToPlot[contactToPlot['part_id'].isin(rsi['part_id'])]
        cts_d = cts.reset_index().groupby('part_id').apply(lambda x: list(x['index']), include_groups=False).to_dict()
        caseid_n = [nnk for nk,k in rsi.set_index('new_id')['part_id'].to_dict().items() for nnk in [nk]*len(cts_d[k])]
        cts = cts.loc[[cv for c in rsi['part_id'] for cv in cts_d[c]]].copy()
        cts['part_id'] = caseid_n
        rsi['old_id'] = rsi['part_id'].values
        rsi['part_id'] = rsi['new_id'].values
        ### Check if all groups are present in future-index
        if group_loc is None:
            group_loc = contactToPlot[[group_col,order_col,group_pop]].dropna().drop_duplicates().sort_values(order_col).loc[:,group_col].dropna().unique()
        cts_index_group = cts[[group_col,order_col,group_pop]].dropna().drop_duplicates().sort_values(order_col).loc[:,group_col].dropna().unique()
        # relative presence of population group_col within sample
        weight = (rsi.groupby(group_col)['part_id'].count()/rsi.shape[0])
        weight = weight.reindex(index=group_loc)

        added_g = []
        for g in group_loc:
            if g in cts_index_group:
                continue
            else:
                added_g.append(g)
                # rsi, cts = add_fake_diag(g, cts, rsi, group_col, group_pop, order_col, rsi_all)

        # compute reciprocal and return avg.contacts per group_col
        sym, tmp1, Q = symmetrize_matrix(cts,rsi,group_col,group_pop, order_col, index_NA=index_NA)
        # compute assortativity index
        tmp1 = tmp1.reindex(index=group_loc, columns=group_loc)
        Q1g = []
        for g in group_loc:
            if not g in added_g:
                try:
                    tmp1.loc[:,g] = tmp1.loc[:,g].fillna(0)
                    t1g = (tmp1.loc[g,g]/tmp1.loc[:,g].sum() - 1/(tmp1.loc[:,g].dropna().shape[0]))    
                    Q1g.append(t1g)
                except: print(tmp1, g, added_g); raise Exception('STOP')
            else:
                Q1g.append(np.nan)
        # print(added_g, tmp1)
        # for g in added_g: 
        #     tmp1[g] = np.nan
                # sort tmp matrix
        tmp1 = tmp1.loc[group_loc[::-1],group_loc]
        if save_bootstraps: 
            if not os.path.exists(bootstrap_dir):
                os.makedirs(bootstrap_dir)
            tmp1.loc[group_loc,group_loc].to_csv(f'{bootstrap_dir}contact_matrix_{group_col}_bootstrap_{i}.csv')
        # setting up assortativity index and averaage contacts (non symmetrized)
        q1g = pd.Series(Q1g,index=group_loc)
        q1gs.append(q1g)
        tmp1s.append(tmp1)
        Qs.append(Q)
        weights.append(weight)
    if minimal_output: return pd.concat(q1gs, axis=1), tmp1s
    else:              return pd.concat(q1gs, axis=1), tmp1s, Qs, weights

#####################################################################
#####################################################################
#####################################################################

####################################################################################
def plot_contact_matrix(tmpx1, q1g, dsx1, fig_title='Contacts', A_ylims=[(-.0,0.75),(-0.0,18.5)], C_lims=(0,None), annotations=True, annotation_fontsize=8, square=False, title=None, Q=True, axs=None, matrix_only=False, cbar=True, cbar_loc='top'):
    """
        Plot contact matrix with
        'A': average Q-index and average contacts
        'C': tmpx1 contact matrix
    """
    check = False
    if axs is None:
        check = True
        fig, ax = plt.subplot_mosaic([['A'],['C']], gridspec_kw={'height_ratios': [0.2,.8]}, figsize=(6,9))
        axs = [ax['A'],ax['C']]
    # elif (len(axs)==1)&(matrix_only):
    #     # fig, ax = plt.subplots(figsize=(8,3))
    #     axs = [axs[0]]

    cmap = sns.color_palette("icefire", 50)[:]
    if not matrix_only:
        #### Top panels
        q_color = cmap[40]
        c_color = cmap[10]
        axAt = axs[0].twinx()
        axs[0].tick_params(axis='y', colors=q_color)
        axAt.tick_params(axis='y', colors=c_color)
        try:
            lp1 = sns.lineplot(data=dsx1.reset_index(), x='respondent', y=0, color=c_color, ax=axAt, n_boot=1,errorbar=('sd', 1.96))
            lp2 = sns.lineplot(data=q1g.reset_index(), x='index', y=0, ls='--', color=q_color, ax=axs[0], n_boot=1,errorbar=('sd', 1.96))
            lp1.grid(False)
            lp2.grid(False)
        except:
            print('No line to plot',dsx1.reset_index(),q1g.reset_index())
            axs[0].grid(False)
        # Bottom heatmap panels
    sns.heatmap(tmpx1, cmap=cmap, vmin=C_lims[0], vmax=C_lims[1], annot=annotations, square=square, fmt='.1f', ax=axs[-1], 
                cbar_kws={'label': 'number of contacts','location':cbar_loc},annot_kws={"size":annotation_fontsize}, cbar=cbar)

    # Set haestetics
    if not matrix_only:
        if Q: axs[0].set_title(f'{fig_title} - Q: {round(q1g.mean(),2)}', fontsize=12, fontstyle='italic')
        else: axs[0].set_title(f'{fig_title}', fontsize=12, fontstyle='italic')
        axs[0].set_xlabel('respondent age-group')
        axs[0].set_ylabel('assortativity', color=q_color)
        axAt.set_ylabel('number of contacts', color=c_color)
        axs[0].set_ylim(A_ylims[0])
        axAt.set_ylim(A_ylims[1])
        axs[0].tick_params(axis='x', labelrotation=90)
        axs[0].grid(alpha=0)
    
    axs[-1].set_xlabel('respondent age-group')
    axs[-1].set_ylabel('contact age-group')

    axs[-1].tick_params(axis='x', labelrotation=90)
    axs[-1].tick_params(axis='y', labelrotation=0)
    if not matrix_only: plt.subplots_adjust(hspace=0.2)
    else: plt.subplots_adjust(hspace=0.35)
    
    if matrix_only: return axs[-1]
    elif check: return fig
    else: return axs


def censor_contacts(contacts, censor=100):
    return contacts.groupby('part_id',as_index=False).apply(lambda x: x.head(censor))

def nan_col(x, loc, col):
    if x[col] != loc: return np.nan
    else:             return x['cont_id']
