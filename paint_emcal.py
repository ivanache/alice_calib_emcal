#!/usr/bin/env python

def set_root_style():
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetOptStat(0)

def multipanel_pad(canvas, pad_layout_unscaled, margin = None):
    pad = []
    for i in range(2):
        s = sum(pad_layout_unscaled[i])
        for j in range(len(pad_layout_unscaled[i])):
            pad_layout[i][j] = pad_layout_unscaled[i][j] / s
    pad_layout_1_reversed = list(reversed(pad_layout[1]))
    if margin == None:
        left_margin = canvas.GetLeftMargin()
        right_margin = canvas.GetRightMargin()
        top_margin = canvas.GetTopMargin()
        bottom_margin = canvas.GetBottomMargin()
    else:
        left_margin, right_margin, top_margin, bottom_margin = magin
    lr = left_margin + right_margin
    tb = top_margin + bottom_margin
    for i in range(len(pad_layout[0])):
        x0 = sum(pad_layout[0][:i]) * (1 - lr)
        x1 = lr + sum(pad_layout[0][:i+1]) * (1 - lr)
        for j in reversed(range(len(pad_layout[1]))):
            y0 = sum(pad_layout_1_reversed[:j]) * (1 - tb)
            y1 = tb + sum(pad_layout_1_reversed[:j+1]) * (1 - tb)
            name = 'pad%d' % len(pad)
            canvas.cd()
            pad.append(ROOT.TPad(name, name, x0, y0, x1, y1))
            pad[-1].SetLeftMargin(
                left_margin / (lr + pad_layout[0][i] * (1 - lr)))
            pad[-1].SetRightMargin(
                right_margin / (lr + pad_layout[0][i] * (1 - lr)))
            pad[-1].SetTopMargin(
                top_margin / (tb + pad_layout_1_reversed[j] * (1 - tb)))
            pad[-1].SetBottomMargin(
                bottom_margin / (tb + pad_layout_1_reversed[j] * (1 - tb)))
            pad[-1].SetFillStyle(0)
            pad[-1].SetFillColor(0)
            pad[-1].Draw()
            pad[-1].cd()
    return pad

def to_sm_ieta_iphi(n):
    if n < 11520:
        sm = n / 1152
        n0 = sm * 1152
        nphi = 24
    elif n < 12288:
        sm = 10 + (n - 11520) / 384
        n0 = 11520 + (sm - 10) * 384
        nphi = 8
    elif n < 16896:
        sm = 12 + (n - 12288) / 768
        n0 = 12288 + (sm - 12) * 768
        nphi = 24
    else:
        sm = 18 + (n - 16896) / 384
        n0 = 16896 + (sm - 18) * 384
        nphi = 8
    n1 = n - n0
    ieta = 2 * (n1 / (2 * nphi)) + 1 - (n1 % 2);
    iphi = (n1 / 2) % nphi;
    return sm, ieta, iphi

def to_eta_phi(sm, ieta, iphi):
    coeff_eta = [
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538955994897965, 0.0012327876984126996,
          -0.00005219165046974383, -0.013757946609494863,
          -0.000025332446808511946, 1.086380600873213e-6 ],
        [ -0.0072721088435373045, -0.00004216269841272758,
          1.131762228686001e-6, -0.01375794660949486,
          -0.00002533244680851135, 1.0863806008729037e-6 ],
        [ 0.6537747727272758, 0.0012291942148758277,
          -0.00005344322673373451, -0.01375054364989851,
          -0.00002391323211245261, 1.0397057440196336e-6 ],
        [ -0.22750791958042052, -0.0004878840193896497,
          0.000021212348669116138, -0.01375054364989848,
          -0.00002391323211246439, 1.0397057440201827e-6 ],
        [ 0.6537747727272758, 0.0012291942148758277,
          -0.00005344322673373451, -0.01375054364989851,
          -0.00002391323211245261, 1.0397057440196336e-6 ],
        [ -0.22750791958042052, -0.0004878840193896497,
          0.000021212348669116138, -0.01375054364989848,
          -0.00002391323211246439, 1.0397057440201827e-6 ],
        [ 0.6537747727272758, 0.0012291942148758277,
          -0.00005344322673373451, -0.01375054364989851,
          -0.00002391323211245261, 1.0397057440196336e-6 ],
        [ -0.22750791958042052, -0.0004878840193896497,
          0.000021212348669116138, -0.01375054364989848,
          -0.00002391323211246439, 1.0397057440201827e-6 ],
        [ 0.6538955994897965, 0.0012327876984126996,
          -0.00005219165046974383, -0.013757946609494863,
          -0.000025332446808511946, 1.086380600873213e-6 ],
        [ -0.0072721088435373045, -0.00004216269841272758,
          1.131762228686001e-6, -0.01375794660949486,
          -0.00002533244680851135, 1.0863806008729037e-6 ] ]

    coeff_phi = [
        [ 1.4153378750000025, 0.0006568888888901503,
          0.013518083333333313, -0.000057045289855166094 ],
        [ 1.4159947638888941, -0.0006568888888867641,
          0.013461038043478247, 0.00005704528985495655 ],
        [ 1.7644118472222259, 0.000646180555557041,
          0.01351747826086953, -0.00005638586956534543 ],
        [ 1.7650580277777765, -0.0006461805555546355,
          0.013461092391304308, 0.00005638586956517725 ],
        [ 2.1134724027777834, 0.0006524861111138665,
          0.013517768115941999, -0.00005654166666683834 ],
        [ 2.1141248888888984, -0.0006524861111081304,
          0.013461226449275318, 0.000056541666666427527 ],
        [ 2.462534333333343, 0.0006581805555594584,
          0.013518074275362234, -0.00005697644927563384 ],
        [ 2.463192513888897, -0.0006581805555524085,
          0.013461097826086869, 0.000056976449275180296 ],
        [ 2.8116099861111175, 0.0006471666666698442,
          0.013517443840579645, -0.00005650181159440986 ],
        [ 2.8122571527777906, -0.0006471666666621282,
          0.013460942028985465, 0.00005650181159390807 ],
        [ -3.12222569444445, 0.0006597222222129356,
          0.013424007936509517, -0.00005902777777630214 ],
        [ -3.121565972222227, -0.0006597222222309393,
          0.013364980158731608, 0.000059027777779336 ],
        [ -1.7262874166666704, 0.0008678541666655771,
          0.013520961956521787, -0.00007551086956511075 ],
        [ -1.7254195625000033, -0.000867854166667594,
          0.013445451086956564, 0.00007551086956536232 ],
        [ -1.3772168333333357, 0.000859541666666056,
          0.013520576086956563, -0.00007492391304337749 ],
        [ -1.3763572916666709, -0.000859541666667937,
          0.013445652173913077, 0.00007492391304362568 ],
        [ -1.0281504166666693, 0.0008688541666657723,
          0.013520394021739162, -0.00007519021739121045 ],
        [ -1.0272815625000018, -0.000868854166667404,
          0.013445203804347865, 0.00007519021739139448 ],
        [ -0.6787590277777789, 0.0006451388888869277,
          0.013420585317460642, -0.00005367063492026231 ],
        [ -0.67811388888889, -0.0006451388888908697,
          0.01336691468254001, 0.000053670634921005446 ] ]

    if abs(ieta - round(ieta) >= 0.25):
        ieta_mod_2 = 0.5
    else:
        ieta_mod_2 = ieta % 2
    eta = coeff_eta[sm][0] + coeff_eta[sm][1] * iphi + \
          coeff_eta[sm][2] * iphi * iphi + \
          coeff_eta[sm][3] * ieta + \
          coeff_eta[sm][4] * ieta * iphi + \
          coeff_eta[sm][5] * ieta * iphi * iphi
    phi = coeff_phi[sm][0] + coeff_phi[sm][1] * ieta_mod_2 + \
          coeff_phi[sm][2] * iphi + \
          coeff_phi[sm][3] * iphi * ieta_mod_2
    if phi < -2:
        phi += 2 * math.pi
    return eta, phi

def update(canvas, pad, root_histogram, root_histogram_index = 0):
    list_index = None
    try:
        index_list = (i for i, v in enumerate(root_histogram)
                      if isinstance(v, list)).next()
    except StopIteration:
        None
    if list_index != None:
        del root_histogram[list_index:]
    root_histogram.append([])
    azimuth_range = [(-1.85, -0.5), (1.3, 3.35)]
    for i in range(2):
        root_histogram[-1].append(ROOT.TH2D(
            'root_histogram%d_%d' %
            (len(root_histogram), len(root_histogram[-1])), '',
            1, azimuth_range[i][0], azimuth_range[i][1],
            1, -0.75, 0.75))
        if i == 0:
            root_histogram[-1][-1].GetXaxis().SetLabelFont(42)
            root_histogram[-1][-1].GetYaxis().SetLabelFont(42)
            root_histogram[-1][-1].GetYaxis().SetTitle('\\eta')
        elif i == 1:
            root_histogram[-1][-1].Fill(0, 0, -1)
            root_histogram[-1][-1].GetXaxis().SetLabelFont(42)
            root_histogram[-1][-1].GetXaxis().SetTitle('\\varphi')
            root_histogram[-1][-1].GetYaxis().SetLabelSize(0)
            root_histogram[-1][-1].GetZaxis().SetLabelFont(42)
        root_histogram[-1][-1].SetMinimum(
            root_histogram[root_histogram_index].GetMinimum())
        root_histogram[-1][-1].SetMaximum(
            root_histogram[root_histogram_index].GetMaximum())
        root_histogram[-1][-1].GetZaxis().SetRangeUser(
            root_histogram[root_histogram_index].GetMinimum(),
            root_histogram[root_histogram_index].GetMaximum())
    root_histogram.append([])
    ncell = 17664
    ncell_emcal = 12288
    pad[1].cd()
    root_histogram[-2][1].Draw('colz')
    pad[1].Update()
    root_histogram[-2][1].GetListOfFunctions().\
        FindObject('palette').SetX2NDC(0.9125)
    for i in range(ncell):
        if i == ncell_emcal:
            pad[0].cd()
            root_histogram[-2][0].Draw('colz')
            pad[0].Update()
        sm, ieta, iphi = to_sm_ieta_iphi(i)
        p = [
            to_eta_phi(sm, ieta - 0.5, iphi - 0.5),
            to_eta_phi(sm, ieta - 0.5, iphi + 0.5),
            to_eta_phi(sm, ieta + 0.5, iphi + 0.5),
            to_eta_phi(sm, ieta + 0.5, iphi - 0.5)
        ]
        x = array.array('d', [p[k][1] for k in range(4)])
        y = array.array('d', [p[k][0] for k in range(4)])
        root_histogram[-1].append(ROOT.TPolyLine(4, x, y))
        root_histogram[-1][-1].SetFillColor(
            root_histogram[-2][1].GetListOfFunctions().
            FindObject('palette').GetValueColor(
                root_histogram[root_histogram_index].
                GetBinContent(i + 1)))
        root_histogram[-1][-1].SetLineColor(ROOT.kBlack);
        root_histogram[-1][-1].SetLineWidth(1);
        root_histogram[-1][-1].Draw('f');
    canvas.Update()

import os, sys
sys.path.append(os.path.join(os.environ['ROOTSYS'], 'lib'))
import math, array, random, ROOT

application_name = ' '.join(sys.argv)
set_root_style()
canvas = ROOT.TCanvas(
    'canvas%d' % 0, application_name, 1280 + 4, 500 + 33)
canvas.SetLeftMargin(0.03125 * 1.5)
canvas.SetRightMargin(0.03125 * 2.5)
canvas.SetTopMargin(0.03125)
pad_layout = ([1.3, 2.0], [1.0])
pad = multipanel_pad(canvas, pad_layout)

ncell = 17664

root_histogram = []
for i in range(1):
    root_histogram.append(ROOT.TH1D(
        'root_histogram%d' % len(root_histogram), '',
        ncell, -0.5, ncell + 0.5))

for i in range(ncell):
    root_histogram[-1].Fill(i, random.uniform(0, 1))

canvas.cd()
for i in range(len(pad)):
    pad[i].Draw()

update(canvas, pad, root_histogram, 0)

ROOT.gApplication.Run()
