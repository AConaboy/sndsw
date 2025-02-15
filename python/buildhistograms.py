import ROOT

title_dict = {"dxL":["#Delta(x-left,x-track) [cm]", "Counts"], "dxR":["#Delta(x-right,x-track) [cm]", "Counts"],
              "dxB":["#Delta(x-barycentre,x-track) [cm]", "Counts"], "dxLvEx":["#Delta(x-left,x-track) [cm]", "x-track [cm]"],
              "dxRvEx":["#Delta(x-right,x-track) [cm]", "x-track [cm]"], "dxBvEx":["#Delta(x-barycentre,x-track) [cm]", "x-track [cm]"],
              "lambda_x":["#lambda_{x} [cm]", "Counts"], "lambda_xvEx":["#lambda_{x} [cm]", "x-track [cm]"], "lambda_xvQDC":["#lambda_{x} [cm]", "QDC [a.u]"],
              "dyB":["#Delta(y-barycentre,y-track) [cm]", "Counts"], "dyBvEy":["#Delta(y-barycentre,y-track) [cm]", "y-track [cm]"],
              "lambda_y":["#lambda_{y} [cm]", "Counts"], 'lambda_yvEy':["#lambda_{y} [cm]", "y-track [cm]"],
              "relQDC":["QDC_{bar}/QDC_{plane}", "Counts"]}

def MakeReconstructMuonHistograms():

    hist_dict = {}

    ####### HISTOGRAMS FOR X ANALYSIS #######

    for x in ['dxL', 'dxR', 'dxB', 'lambda_x']:

        for xm in ('relQDC', 'maxQDC'):
            for plane in range(5):
                
                # Make 1d hist
                histname = f'{x}_{xm}-plane{plane}'
                title = title_dict[f'{x}'][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[f'{x}'][0]+';'+title_dict[f'{x}'][1]
                hist_dict[histname] = ROOT.TH1F(histname, title, 80, -40, 40)

                # Make 2d hist v Ex
                histname = f'{x}_{xm}vEx-plane{plane}'
                title = title_dict[f'{x}vEx'][0].split(' ')[0]+' v '+title_dict[f'{x}vEx'][1].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[f'{x}vEx'][0]+';'+title_dict[f'{x}vEx'][1]+';Counts'
                hist_dict[histname] = ROOT.TH2F(histname, title, 80, -40, 40, 120, -110, 10)

            if x=='lambda_x': continue

            # Make 1d hist
            histname = f'{x}_{xm}-total'
            title = title_dict[x][0].split(' ')[0]+' for full HCAL;'+title_dict[x][0]+';'+title_dict[x][1]
            hist_dict[histname] = ROOT.TH1F(histname, title, 80, -40, 40)

            # Make 2d hist v Ex
            histname = f'{x}_{xm}vEx-total'
            title = title_dict[f'{x}vEx'][0].split(' ')[0]+' v '+title_dict[f'{x}vEx'][1].split(' ')[0]+' for full HCAL;'+title_dict[f'{x}vEx'][0]+';'+title_dict[f'{x}vEx'][1]
            hist_dict[histname] = ROOT.TH2F(histname, title, 80, -40, 40, 120, -110, 10)                          

    for x in ["relQDC"]:
        for plane in range(5):
            histname=f'{x}-plane{plane}'
            title = title_dict[x][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[x][0]+';'+title_dict[x][1]
            hist_dict[histname] = ROOT.TH1F(histname, title, 50, 0, 1)
        histname = f'{x}-total'
        title = title_dict[x][0].split(' ')[0]+' for full HCAL;'+title_dict[x][0]+';'+title_dict[x][1]
        hist_dict[histname] = ROOT.TH1F(histname, title, 50, 0, 1)

    ####### HISTOGRAMS FOR Y ANALYSIS #######

    for y in ['dyB']:
        for plane in range(5):

            # Making 1d hist
            histname = f'{y}-plane{plane}'
            title = title_dict[y][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[y][0]+';'+title_dict[y][1]
            hist_dict[histname] = ROOT.TH1F(histname, title, 60, -30, 30)

            # Making 2d hist
            histname = f'{y}vEy-plane{plane}'
            title = title_dict[f'{y}vEy'][0].split(' ')[0]+' v '+title_dict[f'{y}vEy'][1].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[f'{y}vEy'][0]+';'+title_dict[f'{y}vEy'][1]
            hist_dict[histname] = ROOT.TH2F(histname, title, 60, -30, 30, 80, 0, 80)

        histname = f'{y}-total'
        title = title_dict[y][0].split(' ')[0]+' for full HCAL;'+title_dict[y][0]+';'+title_dict[y][1]
        hist_dict[histname] = ROOT.TH1F(histname, title, 60, -30, 30)

        histname = f'{y}vEy-total'
        title = title_dict[f'{y}vEy'][0].split(' ')[0]+' v '+title_dict[f'{y}vEy'][1].split(' ')[0]+'for full HCAL;'+title_dict[f'{y}vEy'][0]+';'+title_dict[f'{y}vEy'][1]
        hist_dict[histname] = ROOT.TH2F(histname, title, 60, -30, 30, 80, 0, 80)

    for y in ['lambda_y']:
        for plane in range(5):
            histname = f'{y}-plane{plane}'
            title = title_dict[y][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[y][0]+';'+title_dict[y][1]
            hist_dict[histname] = ROOT.TH1F(histname, title, 100, -10, 90)

    for y in ['lambda_yvEy']:
        for plane in range(5):
            histname = f'{y}-plane{plane}'
            title = title_dict[y][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[y][0]+';'+title_dict[y][1]+';Counts'
            hist_dict[histname] = ROOT.TH2F(histname, title, 100, -10, 90, 60, 15, 65)

    return hist_dict

def MakeShowerProfilesHistograms():

    hist_dict = {}

    ####### HISTOGRAMS FOR X ANALYSIS #######

    for x in ['dxL', 'dxR', 'dxB', 'lambda_x']:

        for xm in ('relQDC', 'maxQDC'):
            # for trackStatus in ('noTrack', 'withTrack'):
            for plane in range(5):
                
                # Make 1d hist
                histname = f'{x}_{xm}-plane{plane}'
                title = title_dict[f'{x}'][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[f'{x}'][0]+';'+title_dict[f'{x}'][1]
                hist_dict[histname] = ROOT.TH1F(histname, title, 80, -40, 40)

                # Make 2d hist
                histname = f'{x}_{xm}vEx-plane{plane}'
                title = title_dict[f'{x}vEx'][0].split(' ')[0]+' v '+title_dict[f'{x}vEx'][1].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[f'{x}vEx'][0]+';'+title_dict[f'{x}vEx'][1]+';Counts'
                hist_dict[histname] = ROOT.TH2F(histname, title, 80, -40, 40, 120, -110, 10)

                if x != 'lambda_x':continue
                # Make 2d hist v QDC
                histname = f'{x}_{xm}vQDC-plane{plane}'
                title = title_dict[f'{x}vQDC'][0].split(' ')[0]+' v '+title_dict[f'{x}vQDC'][1].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[f'{x}vQDC'][0]+';'+title_dict[f'{x}vQDC'][1]+';Counts'
                hist_dict[histname] = ROOT.TH2F(histname, title, 80, -40, 40, 120, 0, 120)

            if x=='lambda_x': continue

            # Make 1d hist
            histname = f'{x}_{xm}-total'
            title = title_dict[x][0].split(' ')[0]+' for full HCAL;'+title_dict[x][0]+';'+title_dict[x][1]
            hist_dict[histname] = ROOT.TH1F(histname, title, 80, -40, 40)

            # Make 2d hist
            histname = f'{x}_{xm}vEx-total'
            title = title_dict[f'{x}vEx'][0].split(' ')[0]+' v '+title_dict[f'{x}vEx'][1].split(' ')[0]+' for full HCAL;'+title_dict[f'{x}vEx'][0]+';'+title_dict[f'{x}vEx'][1]
            hist_dict[histname] = ROOT.TH2F(histname, title, 80, -40, 40, 120, -110, 10)       

            # # Make 2d hist v QDC
            # histname = f'{x}_{xm}vQDC-total'
            # title = title_dict[f'{x}vQDC'][0].split(' ')[0]+' v '+title_dict[f'{x}vQDC'][1].split(' ')[0]+' for full HCAL;'+title_dict[f'{x}vQDC'][0]+';'+title_dict[f'{x}vQDC'][1]
            # hist_dict[histname] = ROOT.TH2F(histname, title, 80, -40, 40, 120, 0, 120)              

    for x in ["relQDC"]:
        for plane in range(5):
            histname=f'{x}-plane{plane}'
            title = title_dict[x][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[x][0]+';'+title_dict[x][1]
            hist_dict[histname] = ROOT.TH1F(histname, title, 51, 0, 1.02)
        histname = f'{x}-total'
        title = title_dict[x][0].split(' ')[0]+' for full HCAL;'+title_dict[x][0]+';'+title_dict[x][1]
        hist_dict[histname] = ROOT.TH1F(histname, title, 50, 0, 1)

    ####### HISTOGRAMS FOR Y ANALYSIS #######

    for y in ['dyB']:
        for plane in range(5):

            # Making 1d hist
            histname = f'{y}-plane{plane}'
            title = title_dict[y][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[y][0]+';'+title_dict[y][1]
            hist_dict[histname] = ROOT.TH1F(histname, title, 60, -30, 30)

            # Making 2d hist
            histname = f'{y}vEy-plane{plane}'
            title = title_dict[f'{y}vEy'][0].split(' ')[0]+' v '+title_dict[f'{y}vEy'][1].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[f'{y}vEy'][0]+';'+title_dict[f'{y}vEy'][1]
            hist_dict[histname] = ROOT.TH2F(histname, title, 60, -30, 30, 80, 0, 80)

        histname = f'{y}-total'
        title = title_dict[y][0].split(' ')[0]+' for full HCAL;'+title_dict[y][0]+';'+title_dict[y][1]
        hist_dict[histname] = ROOT.TH1F(histname, title, 60, -30, 30)

        histname = f'{y}vEy-total'
        title = title_dict[f'{y}vEy'][0].split(' ')[0]+' v '+title_dict[f'{y}vEy'][1].split(' ')[0]+'for full HCAL;'+title_dict[f'{y}vEy'][0]+';'+title_dict[f'{y}vEy'][1]
        hist_dict[histname] = ROOT.TH2F(histname, title, 60, -30, 30, 80, 0, 80)

    for y in ['lambda_y']:
        for plane in range(5):
            histname = f'{y}-plane{plane}'
            title = title_dict[y][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[y][0]+';'+title_dict[y][1]
            hist_dict[histname] = ROOT.TH1F(histname, title, 100, -10, 90)

    for y in ['lambda_yvEy']:
        for plane in range(5):
            histname = f'{y}-plane{plane}'
            title = title_dict[y][0].split(' ')[0]+' for plane '+str(plane+1)+';'+title_dict[y][0]+';'+title_dict[y][1]+';Counts'
            hist_dict[histname] = ROOT.TH2F(histname, title, 100, -10, 90, 60, 15, 65)

    return hist_dict
