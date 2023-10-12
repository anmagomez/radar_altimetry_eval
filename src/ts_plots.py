#Adding the plots to the file 

import matplotlib.dates as mdates
from matplotlib.ticker import StrMethodFormatter
import matplotlib.pyplot as plt
import numpy as np
#Plots
#Plot deviation of the mean from satellite, groound observations and interpolated points
def plot_dev_mean(nfig,dfs, dfg, dfi, source, st_id, altis_name, date_fd, height_fd, labels, labelg, labeli, error_m_text=None, output_f=None):
    plt.figure(nfig, figsize=(20, 6))
    splot=plt.plot(dfs[date_fd], 
                   dfs[height_fd]-np.nanmean(dfs[height_fd]), linestyle='None',marker='o', 
                   markerfacecolor='#20641c',markeredgecolor='#20641c',markersize=5, label=labels)
    gplot=plt.plot(dfg[date_fd], 
                   dfg[height_fd]-dfg[height_fd].mean(), linestyle='solid',
                   markerfacecolor='#0000FF',markeredgecolor='#0000FF',color='#0000FF',marker='o', markersize=2, label=labelg)
    iplot=plt.plot(dfi[date_fd].to_numpy(), 
                   dfi[height_fd]-np.nanmean(dfi[height_fd]), linestyle='None',
                   markerfacecolor='#d95f02', markeredgecolor='#d95f02',marker='s', 
                   markersize=5, label=labeli)
    # plt.plot(altidy, altiwelev, '-ro', markersize=2.5)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    ax=plt.gca()
    if error_m_text is not None:
        plt.text(0.01, 0.8, 'Metrics'+error_m_text, fontsize = 15, transform=ax.transAxes)
        
    plt.xlabel('Time (decimal year)', size=15, weight='bold')
    plt.ylabel('Water elev. dev from mean (m)', size=15, weight='bold')
    plt.grid('on')
    plt.legend(loc='upper right')
    plt.title('Deviation from the mean time series comparison', size=15,
              weight='bold')


    if output_f is not None:
        if not os.path.isdir(output_f):
            os.mkdir(output_f) 
        plt.savefig(output_f+source+'_'+st_id+'_vs_'+altis_name+'.png',
                bbox_inches='tight')
    
    plt.show()
    return plt


#Plots
#Plot deviation of the mean from satellite, groound observations and interpolated points
def plot_dev_mean_publish(fig,nfig,dfs, dfg, source, st_id, altis_name, date_fd, height_fd_g,  height_fd_s, labels, labelg, fill_between=None, error_m_text=None, output_f=None):
    #plt.figure(nfig, figsize=(20, 6))
    
    if 'Sentinel-3A' in altis_name:
        markerfacecolor='#D81B60'
        markeredgecolor='#D81B60'
        marker='o'
        markersize=5
    elif 'Sentinel-3B' in altis_name:
        markerfacecolor='#000000'
        markeredgecolor='#000000'
        marker='o'
        markersize=5
    else:
        markerfacecolor='#004D40'
        markeredgecolor='#004D40'
        marker='*'
        markersize=8

    #Azul#0000FF
    ax=plt.gca()
    ax.axhline(linewidth=1, color='#575656',linestyle='--')
    
    dfg=dfg.sort_values(by=date_fd)

    #Get the time difference in minutes between each data point and extract the median
    temp_diff=dfg[date_fd].diff().apply(lambda x: x/np.timedelta64(1, 'm')).fillna(0).astype('int64')
    mean_diff_time=temp_diff.median()

    #Get the indices in which difference in time between data point is greater than 10 days (14400 min), repeat cycle for Jason 3
    indices = np.where(np.abs(temp_diff.to_numpy()) >= 14400)[0] 

    #Plot only the lines between the intervals of continues measurements
    
    for n, i in enumerate(indices):
        if n == 0:
            plt.plot(dfg[date_fd].to_numpy()[:i], (dfg[height_fd_g]-dfg[height_fd_g].mean()).to_numpy()[:i], linestyle='solid',
                    markerfacecolor='#7570b3',markeredgecolor='#7570b3',color='#7570b3',marker='o', markersize=5)#, label=labelg)
        else:
            plt.plot(dfg[date_fd].to_numpy()[indices[n - 1]:i], (dfg[height_fd_g]-dfg[height_fd_g].mean()).to_numpy()[indices[n - 1]:i], linestyle='solid', markerfacecolor='#7570b3',markeredgecolor='#7570b3',color='#7570b3',marker='o', markersize=5) #, label=labelg)

    #Plot the last part of the series. If not gaps it plots the complete series
    if len(indices)==0:
        i=0
    
    l_idx=len(temp_diff.to_numpy())
    plt.plot(dfg[date_fd].to_numpy()[i:l_idx], (dfg[height_fd_g]-dfg[height_fd_g].mean()).to_numpy()[i:l_idx], linestyle='solid',
                    markerfacecolor='#7570b3',markeredgecolor='#7570b3',color='#7570b3',marker='o', markersize=3) #, label=labelg)
    
    splot=plt.plot(dfs[date_fd], 
                   dfs[height_fd_s]-np.nanmean(dfs[height_fd_s]), linestyle='None',marker=marker, 
                   markerfacecolor=markerfacecolor,markeredgecolor=markerfacecolor,markersize=markersize, label=labels)
    
    # if error_m_text is not None:
    #     plt.text(0.01, 0.8, 'Metrics'+error_m_text, fontsize = 8, transform=ax.transAxes)
        
    # plt.xlabel('Time (decimal year)', size=15, weight='bold')
    # plt.ylabel('Water elev. dev from mean (m)', size=15, weight='bold')
    # plt.grid('on')
    plt.grid(color='#808080', linestyle='-.', linewidth=0.5)
    ymin,ymax=plt.ylim()
    if abs(ymin)>abs(ymax):
        ymax=abs(ymin)+0.1
    elif abs(ymin)<abs(ymax):
        ymin=-abs(ymax)-0.1
        
    #If needed fill between 
    if fill_between is not None:
        fb_fill_fd=fill_between['fill']
        if dfg.loc[dfg[fb_fill_fd].notna()].shape[0]!=0:
            lim=fill_between['lim'] #TODO: Generalize this more
            ax.fill_between(dfg[date_fd], ymin, ymax, where=(dfg[fb_fill_fd] > lim), alpha=0.5, color='#C5C5C5')
        
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 1 decimal places    
    plt.ylim(ymin,ymax)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    for label in ax.get_xticklabels(which='major'):
        label.set(rotation=45, horizontalalignment='right', fontsize=18)
    for ylabel in ax.get_yticklabels(which='major'):
        ylabel.set(fontsize=18)
        
    fig.tight_layout(h_pad=3, w_pad=2)
    # plt.legend(loc='lower right')
    plt.title(labelg+error_m_text, size=20,
              weight='bold')


    if output_f is not None:
        if not os.path.isdir(output_f):
            os.mkdir(output_f) 
        plt.savefig(output_f+source+'_'+st_id+'_vs_'+altis_name+'.png',
                bbox_inches='tight')
    
    return plt