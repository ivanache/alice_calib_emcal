#! /usr/bin/env python
from os import path
from ROOT import TCanvas, TFile, TH1D, TH2D, TLine
from ROOT import gROOT, gSystem, Double
import sys

import bad_channel_lhc15o_3
import bad_channel_lhc15o_mine3
bad_all = sorted(bad_channel_lhc15o_3.bad_all + bad_channel_lhc15o_mine3.hot)

#make some histogram
h2 = TH1D('h2', 'Coeff plot for list 2', 400, -0.5, 3.5)
h3 = TH1D('h3', 'Coeff plot for list 3', 400, -0.5, 3.5)
rat = TH1D('rat', 'Ratio of coeffs between lists 2 and 3', 300, -1.5, 1.5) 
corrG = TH2D('corrG', 'Correlation between the two lists: good channels only', 400, -0.5, 3.5, 400, -0.5, 3.5)
corrG.GetXaxis().SetTitle('c2')
corrG.GetYaxis().SetTitle('c3')
corrB = TH2D('corrB', 'Correlation between the two lists: bad channels only', 400, -0.5, 3.5, 400, -0.5, 3.5)
list2_dict = {}
list3_dict = {}
with open('lhc15o_list2_chiBeforeAfter.txt','r') as list2_txt:
    for line2 in list2_txt:
        cell_id2, eb2, yb2, coeff2, chiA2, chiB2 = [float(i) for i in line2.split()]
        cell_id2 = int(cell_id2)
        eb2 = int(eb2)
        yb2 = int(yb2)
        list2 = [eb2, yb2, coeff2, chiA2, chiB2]
        list2_dict[cell_id2] = list2
        h2.Fill(coeff2)
        #print cell_id2, list2
        #if cell_id2 == 20:
        #    break

#print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

with open('lhc15o_list3_chiBeforeAfter.txt', 'r') as list3_txt:
    for line3 in list3_txt:      
        cell_id3, eb3, yb3, coeff3, chiA3, chiB3 = [float(j) for j in line3.split()]
        cell_id3 = int(cell_id3)
        eb3 = int(eb3)
        yb3 = int(yb3)
        list3 = [eb3, yb3, coeff3, chiA3, chiB3]
        list3_dict[cell_id3] = list3
        h3.Fill(coeff3)
        #print cell_id3, list3
        #if cell_id3 == 20:
        #    break

bad_channels = []
missing_channels = []
outliers = []
for key in list2_dict.keys():
    if key in list3_dict:
        l2 = list2_dict[key]
        l3 = list3_dict[key]
        c2 = float(l2[2])
        c3 = float(l3[2])
        if key in bad_all:
            corrB.Fill(c2, c3)
            bad_channels.append(key)
            #print key, 'is a bad channel'
        else:
            if (c2 < 0.85 or c2 > 1.25) and (1.3 > c3 and c3 > 0.8):
                outliers.append([key, c2, c3])
            if l2[3] < 10:
                corrG.Fill(c2, c3)
                rat.Fill(c2-c3)
            #if (c2
            
    else:
        #print key, 'is missing list3_dict'
        missing_channels.append(key)
        
lCorr = TLine(0,0,3.5,3.5)
#Creating a Canvas and Drawing
cg = TCanvas('cg', 'cg', 600, 400)
cg.SetLogz()
corrG.Draw('colz')
lnG = TLine(0,1,3.5,1)
lnG.Draw()
lCorr.Draw()

cb = TCanvas('cb', 'cb', 600, 400)
cb.SetLogz()
corrB.Draw('colz')
lnB = TLine(1,0,1,3)
lnB.Draw()
lCorr.Draw()

cr = TCanvas('cr', 'cr', 600, 400)
rat.Draw()
#rat.Fit('gaus','','',0.85,1.15)
rat.GetXaxis().SetTitle('c2-c33')



#print bad_channels, '\n\n'
print missing_channels

for n in outliers:
    print n[0]

# ---------------------------------------------------------------------------------- 
raw_input('Press Enter to exit')
