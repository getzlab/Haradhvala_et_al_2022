import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from matplotlib.patches import Patch


def volcano(df,thresh=.25,ax=None,label=False,min_size=1,max_size=50,palette=None):
    gene_show = df[df['adj.P.Val']<thresh].sort_values('logFC').index
    
    df['sig']='nonsignificant'
    df.loc[(df['adj.P.Val']<thresh)&(df['logFC']>0),'sig']='up'
    df.loc[(df['adj.P.Val']<thresh)&(df['logFC']<0),'sig']='down'
    
    if palette is None:
        palette = {'nonsignificant':'gray','up':'blue','down':'red'}
    
    
    df['nlogp'] = -1*np.log10(df['P.Value'])
    df['size'] = df['nlogp']*df['logFC'].abs()
    df['size'] = df['size']/df['size'].max()
    
    if ax is None:
        f,ax = plt.subplots(1)
    g=sns.scatterplot(x="logFC",y="nlogp",data=df,hue='sig',size="size",sizes=(min_size,max_size),
                      ax=ax,palette=palette,legend=False,rasterized=True)
    
    
    if label:
        texts = [ax.text(row['logFC'],row['nlogp'],gene) for 
                 gene,row in df[df['adj.P.Val']<thresh].iterrows()]
        adjust_text(texts)
    ax.grid(False)
    
    
def sankeyish_plot(node_data,edge_data,colors,hatches=None,ax=None,node_height=1,
                   node_width=.1,linewidth=1,divider_color='k',boxes = [],box_width=3,edge_alpha=1):
    
    # node data is list of dictionaries, which each has keys
        # xpos
        # ypos
        # freqs (list of frequencies)
        
    if ax is None:
        f,ax=plt.subplots(figsize=(10,5))
        
    if hatches is None:
        hatches = [None for c in colors]
    
    ### Plot bar of frequencies for each node
    for node in node_data:
        # Determine the bottom position of each subtype
        node['bottom'] = np.hstack([0,np.cumsum(node['freqs'])[0:-1]]) + node['ypos'] - node_height/2

        for bottom,freq,color,hatch in zip(node['bottom'],node['freqs'],colors,hatches):
            ax.bar(node['xpos'],
                   bottom=bottom,
                   height=freq,
                   color=color,
                   width=node_width,edgecolor=divider_color,
                   linewidth=linewidth,
                  hatch=hatch)
        
        # Put boxes around heighlighted populations
        for box in boxes:
            bstart,bend = box
            box_bottom = node['bottom'][bstart]
            box_top = node['bottom'][bend]+node['freqs'][bend]
            box_left = node['xpos']-node_width/2
            box_right = node['xpos']+node_width/2
            
            for xpos in [box_left,box_right]:
                ax.plot([xpos,xpos],[box_bottom,box_top],color='k',linewidth=box_width)
            for ypos in [box_bottom,box_top]:
                ax.plot([box_left,box_right],[ypos,ypos],color='k',linewidth=box_width)
            
    # Plot edges
    for edge in edge_data:
        n1 = node_data[edge[0]]
        n2 = node_data[edge[1]]
        xpos = [n1['xpos']+(node_width/2)*.95,
                    n2['xpos']-node_width/2]
        
        for b1,b2,f1,f2,color,hatch in zip(n1['bottom'],n2['bottom'],n1['freqs'],n2['freqs'],colors,hatches):
            ax.fill_between(xpos,[b1,b2],
                            [b1+f1,b2+f2],color=color,
                            edgecolor=divider_color,
                            linewidth=linewidth,alpha=edge_alpha,hatch=hatch)
    return(ax)
        
    
    
    
def evolution_plot(F,ax=None,legend=True,**kwargs):
    
    if ax is None:
        f,ax=plt.subplots(figsize=(10,5))
    
    tps = ['Baseline','Infusion','D7','D7-CAR-T']
    xpos = [0,.5,1,1]
    ypos = [0,-1,1,-1]

    cd4_hatch = '/'
    ncd4 = F.loc['CD4 T'].shape[0]

    spacer_width=.1

    class_colors = {'CM':'purple','EM':'orange','EMRA':'sienna',
                    'Naive':'gray','T-reg':'red','CD4+ CTL':'teal'}

    colors = [class_colors[cl] for (st,cl) in F.index]

    # Create plot data
    node_data = [{'xpos':xpos[i],'ypos':ypos[i],'freqs':F[tps[i]]} for i in range(0,len(tps))]
    
    colors.insert(ncd4,'white')

    node_data = [{'xpos':xpos[i],'ypos':ypos[i],'freqs':np.insert(F[tps[i]].values,ncd4,spacer_width)} \
             for i in range(0,len(tps))]

    edge_data = [(1,3),(0,2),(0,1)]

    hatches = [cd4_hatch if ind[0]=='CD4 T' else None for ind in F.index]

    # Make plot
    sankeyish_plot(node_data,edge_data,colors,hatches=hatches,ax=ax,**kwargs)

    # Format axes
    ax.set_xticks([0,.5,1]);
    ax.set_xticklabels(['Baseline','Infusion','Day 7'],fontsize=16)
    ax.set_yticks([])
    ax.text(-.08,.3,'CD8 T',rotation=90,va='center',ha='center',fontsize=14)
    ax.text(-.08,-.3,'CD4 T',rotation=90,va='center',ha='center',fontsize=14)

    ax.text(1.08,-1,'CAR+',rotation=-90,va='center',ha='center',fontsize=14)
    ax.text(1.08,1,'CAR-',rotation=-90,va='center',ha='center',fontsize=14)
    
    if legend:
        # Make legend
        legend_size=16

        hatch_dict = {'CD4 T':cd4_hatch + cd4_hatch,'CD8 T':None}

        legend_elements = [Patch(facecolor=class_colors[cl],
                         label=cl) for cl in class_colors.keys()]
        cd48_elements = [Patch(facecolor='grey',
                         label=subtype,hatch=hatch_dict[subtype]) for subtype in ['CD4 T','CD8 T']]

        # Create and add subtype legend
        first_legend = ax.legend(handles=legend_elements, loc='lower left',
                             prop={'size': legend_size},bbox_to_anchor=(1.01,0))
        ax.add_artist(first_legend)

        # Create CD4/8 legend
        ax.legend(handles=cd48_elements,prop={'size': legend_size}, loc='upper left',bbox_to_anchor=(1.01,1))



def plot_bracket(x,c,t=None,direction=1,root_x=0,root_y=0,
                 yspan=.5,xspan=1,size_scale=1000,size_min=0,ax=None,fontsize=12):
    
    if ax is None:
        f,ax = plt.subplots(1)
    
    xshift = xspan*direction
    
    # Bracket base 
    ax.plot([root_x,root_x+xshift/2],[root_y,root_y],color='k')
    
    # Place brack ends evenly spaced
    yspace=np.linspace(root_y-yspan/2,root_y+yspan/2,len(x))
    
    # Vertical bar
    ax.plot([root_x+xshift/2,root_x+xshift/2],[yspace[0],yspace[-1]],color='k')
    
    for i in range(0,len(x)):
        
        
        
        if type(x[i]) is tuple:
            # Recursively plot other brackets
            plot_bracket(x[i][1],x[i][2],x[i][3],direction,root_x+xshift,yspace[i],
                         yspan=yspan/len(x),xspan=xspan,size_scale=size_scale,ax=ax,fontsize=fontsize)
            xval=x[i][0]
        else:
            xval=x[i]
        

        
        # Text
        if not t is None:
            ax.text(root_x+xshift*2,yspace[i],t[i],color=c[i],ha='center',va='center',fontsize=fontsize)
            
        # Forks of bracket
        ax.plot([root_x+xshift/2,root_x+xshift],[yspace[i],yspace[i]],color='k')
        
        # Bubble at end of bracket
        ax.scatter(root_x+xshift,yspace[i] ,s=xval*size_scale+size_min,edgecolors='k',c=c[i], zorder=2.5)
            
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
            
def pat_breakdown_plot(X,key='name',direction=1,
                       xspan=1,yspan=.5,size_scale=1000,
                       root=(0,0),norm_total=True,ax=None,label=False,fontsize=12):
    
    
    if ax is None:
        f,ax=plt.subplots(1)

    
    subtypes = ['CD4 T','CD8 T']
    subtype_colors = {'CD4 T':'blue','CD8 T':'orange'}
    class_colors = {'T-reg':'green','EM':'tab:olive','SLEC':'grey','CM':'tab:purple','CTL':'teal'}
    
    
    # Accumulate data to plot
    data = list()
    
    for subtype,g in X.groupby('subtype'):
        
        nsubtype = g.groupby('barcode')['n'].sum()
        ntot = X.groupby('barcode')['n'].sum()

        
        nsubtype[~nsubtype.index.isin(ntot.index)]=0
        
        subtype_mean,sd = est_frequency_center_and_sd(nsubtype,ntot.loc[nsubtype.index],est_type='median')
        

        Xf = pd.Series({clust:est_frequency_center_and_sd(gs['n'].values,
                                ntot.loc[gs['barcode']].values,est_type='median')[0] for clust,gs in g.groupby(key)})
        # For each subcluster
        
        Xf = Xf.iloc[::-1]
        
        class_map = g[[key,'class']].drop_duplicates().set_index(key)['class']
        
        clust_colors = [class_colors[class_map[clust]] for clust in Xf.index]
        
        
        # Add for later plotting
        #labels=[class_map[i] for i in Xf.index] if label else None
        labels=[i for i in Xf.index] if label else None

        data.append((subtype_mean,Xf,clust_colors,labels))
        
    plot_bracket(data,["blue","orange"],direction=direction,root_x=root[0],root_y=root[1],
                 xspan=xspan,yspan=yspan,size_scale=size_scale,ax=ax,fontsize=fontsize)
                    
                    
                    
from scipy.stats import betabinom
from scipy.stats import beta
from scipy.optimize import minimize

def bbloss(X,N,a,b):
    return(-1*sum(betabinom.logpmf(X,N,a,b)))

def est_frequency_center_and_sd(x,N,est_type='mean'):
    # Fit a beta binomial to data
    fit = minimize(lambda ab: bbloss(x,N,ab[0],ab[1]),(1,1),bounds=[(1,None),(1,None)])
    a,b = fit.x

    if est_type=='mean':
        center = a/(a+b)
    elif est_type=='median':
        center = beta.ppf(.5,a,b)
    else:
        raise("Error: undefined central tendancy")
        
    sd = np.sqrt(a*b/((a+b)**2*(a+b+1)))
    
    return(center,sd)

def make_circle_legend(entries,scale,ax,**kwargs):
    
    gs = list()
    for entry in entries:
        g = ax.scatter([],[], s=scale*entry, marker='o', color='#555555')
        gs.append(g)

    ax.legend(gs,
       [f'{entry*100}%' for entry in entries],
       scatterpoints=1,
       loc='upper right',
       fontsize=8,frameon=False,**kwargs)

def clone_tracking_plot(adatas,product,subtype,key='name',ax=None,max_clones=10,sscale=500):

    class_colors = {'T-reg':'green','EM':'tab:olive','SLEC':'grey','CM':'tab:purple','CTL':'yellow'}
    response_palette = {'R':'teal','NR':[1,.2,0]}
 
    if ax is None:
        f,ax = plt.subplots(1,figsize=(7,7))

    tps = ['Infusion','D7-CAR-T']
    
    tp_pos = {tps[0]:0,tps[1]:1}
    
    adict = {tp: adatas[product][subtype][tp] for tp in tps}
    
    # Count tcrs for each patient+tcr  x cluster
    Xdict = dict()
    for tp in tps:
        idx = ~adict[tp].obs['tcr'].str.match('.*nan')
        Xdict[tp] = pd.crosstab([adict[tp].obs.loc[idx,'tcr'],
                                 adict[tp].obs.loc[idx,'barcode']],
                                                          adict[tp].obs.loc[idx,key]).iloc[:,::-1]
        
    # Find most frequent TCRs
    tcr_sort = Xdict[tps[-1]].sum(axis=1).sort_values(ascending=False).reset_index().rename(columns={0:'nd7'})
    
    # Number of subtypes for each timepoint
    ns = {tp:len(Xdict[tp].columns) for tp in tps}


    # Bottom of y-range for each subtype
    class_pos = {tp:dict(zip(Xdict[tp].columns,np.linspace(0,1,ns[tp],endpoint=False))) for tp in tps}
    
    # Phase x cluster count matrix
    phase_counts = {tp:pd.crosstab(adict[tp].obs['phase'],adict[tp].obs[key]).iloc[:,::-1] for tp in tps}


    tcr_totals = {tp:Xdict[tp].reset_index().groupby('barcode').sum().sum(axis=1) for tp in tps}
    xcurve = np.linspace(1e-3,1,endpoint=False)
    
    tdfs = dict()
    for tp in tps:
        tcr_sort = pd.DataFrame(Xdict[tp].sum(axis=1).sort_values(ascending=False),columns=['n'])
        tcr_sort = tcr_sort.reset_index()
        tcr_sort['f'] = tcr_sort['n']/[tcr_totals[tp][b] for b in tcr_sort['barcode']]
        tcr_sort = tcr_sort.set_index(['tcr','barcode'])

        tcr_sort['class'] = Xdict[tp].loc[tcr_sort.index].idxmax(axis=1)
    
        tdfs[tp] = tcr_sort

    tdf = tdfs[tps[0]].join(tdfs[tps[1]],lsuffix='_first',rsuffix='_second',how='inner')  
    
    tdf = tdf.sort_values('n_second',ascending=False)

    tdf = tdf.head(max_clones)
    
    tdf['rank_first'] = tdf['f_first'].rank(method='first')
    tdf['rank_second'] = tdf['f_second'].rank(method='first')
    
    # Set up destination lines and cell cycle pie charts
    for tp in tps:
        # Number of subtypes
        n = ns[tp]
        for i in range(0,ns[tp]):
            ax.plot([tp_pos[tp],tp_pos[tp]],[i/n,i/n+.9/n],zorder=0,color='gray',linewidth=3)
            ha = 'right' if tp==tps[0] else 'left'
            ax.text(-.2 + 1.4*tp_pos[tp],(i+.5)/ns[tp],
                    Xdict[tp].columns[i],
                    ha=ha,va='center')
            
            ax.pie(phase_counts[tp].iloc[:,i],center=[-.1 + 1.2*tp_pos[tp],(i+.5)/n],
                colors = [[0,0,.3],[1,0,.2],[1,.2,.2]],radius=.3/4)

    # Plot clones
    for (tcr,pat),row in tdf.iterrows():
        if '-R-' in pat:
            color=response_palette['R']
        else:
            color=response_palette['NR']
            
        ypos_first = class_pos[tps[0]][row['class_first']] + row['rank_first']/tdf.shape[0]/ns[tps[0]] * .9
        ypos_second = class_pos[tps[-1]][row['class_second']] + row['rank_second']/tdf.shape[0]/ns[tps[-1]] * .9
        
        ax.scatter([0],[ypos_first],s=row['f_first']*sscale,color=color,edgecolors='k',linewidth=.5)
        ax.scatter([1],[ypos_second],s=row['f_second']*sscale,color=color,edgecolors='k',linewidth=.5)
        ax.plot(xcurve,ypos_first + (ypos_second-ypos_first)/(1+(xcurve/(1-xcurve))**-3),color=color)
    

    
    # Make legend
    gll = ax.scatter([],[], s=sscale*.25, marker='o', color='#555555')
    gl = ax.scatter([],[], s=sscale*.1, marker='o', color='#555555')
    ga = ax.scatter([],[], s=sscale*.01, marker='o', color='#555555')

    ax.legend((gll,gl,ga),
       ('25%', '5%', '1%'),
       scatterpoints=1,
       loc='upper right',
       fontsize=8,frameon=False)

    #ax.grid(False)
    ax.axis('off')
    ax.set_ylim(-.1,1.1)


