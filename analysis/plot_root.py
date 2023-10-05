from sys import platform
from typing import Optional
import numpy as np
import pandas as pd
if platform == "linux" or platform == "linux2":
    import ROOT
    # %jsroot off
    
def plot_hist(data:pd.DataFrame, x, labels:Optional[dict]=None, colorpalette=None, bins=128, xlabel:str="", ylabel:str="", xunits:str="", yunits:str="", normalize=False, title="Some Plot"):
    """Plots one ore multiple histograms using ROOT

    Args:
        data (pd.DataFrame): dataframe, i.e. dict with of structure: keys als column names and lists as values
        x (_type_): columns
        labels (_type_, optional): _description_. Defaults to None.
        colorpalette (_type_, optional): If None, defaults to a predefined array of colors. See https://root.cern.ch/doc/master/classTColor.html#C01. Defaults to None.
        bins (int, optional): number of bins. Defaults to 128.
        xlabel (str, optional): label for x axis. Defaults to "".
        ylabel (str, optional): label for y axis. Defaults to "".
        xunits (str, optional): units for x axis. Defaults to "".
        yunits (str, optional): units for y axis. Defaults to "".
        normalize (bool, optional): whether or not to normalize by function area. Defaults to False.
        title (str, optional): _description_. Defaults to "Some Plot".

    Returns:
        _type_: _description_
    """
    
    # Get colorpalette
    pal = colorpalette if colorpalette is not None else [2,4,5,6,7,8,9,10,11]

    # Autoscaling
    g_min = 0.98*data[x].min().min()
    g_max = 1.02*data[x].max().max()
    
    canv = ROOT.TCanvas("c_name", "c_title", 800, 600)
    
    # Create histograms
    hists = []
    for i in range(0, len(x)):
        column = x[i]
        h_title = title if i == 0 else ""
        h_name  = labels[i] if labels is not None else column
        
        hists.append(ROOT.TH1D(h_name, h_title, bins, g_min, g_max))

    xaxis = hists[0].GetXaxis()
    xaxis.SetTitle(xlabel + (" [" + xunits + "]" if xunits != "" else ""))
    
    ylabel = "Normalized" if (normalize and ylabel == "") else ""
    yaxis = hists[0].GetYaxis()
    yaxis.SetTitle(ylabel + (" [" + yunits + "]" if yunits != "" else ""))
    
    #ROOT.gStyle.SetErrorX(0)

    # Assume hist_zhh and hist_zzh are of same size, but skip NaN values
    for j in range(0, len(x)):
        non_nan = data[x[j]]
        non_nan = non_nan[~np.isnan(non_nan)]
        for i in range(0, len(non_nan)):
            hists[j].Fill(non_nan[i])
    
    # Get maximum bin value for scaling
    maxval = 0
    if normalize:
        for hist in hists:
            maxval = max(maxval, hist.GetMaximum()/hist.Integral())
            hist.Scale(1/hist.Integral())
    else:
        for hist in hists:
            maxval = max(maxval, hist.GetMaximum())
        
    hists[0].SetAxisRange(0., 1.02*maxval, "Y")

    # Adjust histogram styling
    for i in range(0, len(hists)):        
        hist  = hists[i]
        color = pal[i]
        
        hist.SetFillStyle(3004)
        hist.SetFillColorAlpha(color, 0.35)
        hist.SetLineWidth(2)
        hist.SetLineColor(color)
    
    # Draw histograms
    for i in range(0, len(hists)):        
        hist  = hists[i]
        hist.Draw("HIST" if i == 0 else "HISTSAMES")
    
    canv.Draw()
    
    ROOT.gPad.Update()

    # Adjust styling of statsboxes
    for i in range(0, len(hists)):
        hist = hists[i]
        sb = hist.FindObject("stats")
        
        hist.GetListOfFunctions().Remove(sb)
        hist.SetStats(0)
    
        sb.SetX1NDC(0.7)
        sb.SetX2NDC(0.95)
        sb.SetY1NDC(0.75 - i*0.15)
        sb.SetY2NDC(0.88 - i*0.15)
        sb.SetLineColor(pal[i])
        sb.SetTextColor(pal[i])
    
        sb.Draw()

    canv.Modified()
    canv.Update()
    
    root_ref = [canv, hists]

    return root_ref