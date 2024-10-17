import ROOT

min_clock_cycle = -0.5
max_clock_cycle = 1.2

scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
muFilterDet = ROOT.gROOT.GetListOfGlobals().FindObject('MuFilter')

a = ROOT.TVector3()
b = ROOT.TVector3()

def selectHits(event, MC=False, return_tref=False):
#    return event.Digi_ScifiHits

    h = []
    for i_orientation in range(2):
        for i_station in range(5):
            h.append(ROOT.TH1D("hHitTime_{}_{}".format(i_orientation, i_station), ";SciFi hit time [ns]", 200, 0., 100.))

    if MC:
        TDC2ns = 1
    else:
        TDC2ns = 1E9/160.316E6

    # Try station by station
    for hit in event.Digi_ScifiHits:
        if not hit.isValid():
            continue
        time = hit.GetTime()*TDC2ns
        if not MC:
            time = scifiDet.GetCorrectedTime(hit.GetDetectorID(), time,0)
        this_orientation = 0
        if hit.isVertical():
            this_orientation = 1
        this_station = hit.GetStation() - 1
        h[this_orientation*5+this_station].Fill(time)

    peakTime = []
    for i_orientation in range(2):
        for i_station in range(5):
            i = i_orientation*5 + i_station
            peakTime.append((h[i].GetMaximumBin()-0.5)*(h[i].GetXaxis().GetXmax()-h[i].GetXaxis().GetXmin())/h[i].GetNbinsX() + h[i].GetXaxis().GetXmin())
    
    sel_hits = []
    for hit in event.Digi_ScifiHits:
        if not hit.isValid(): 
            continue
        time = hit.GetTime()*TDC2ns
        if not MC:
            time = scifiDet.GetCorrectedTime(hit.GetDetectorID(), time,0)

        this_orientation = 0
        if hit.isVertical():
            this_orientation = 1
        this_station = hit.GetStation()-1
        this_i = this_orientation*5+this_station
        
        if time > peakTime[this_i] + max_clock_cycle*1E9/160.316E6:
            continue
        if time < peakTime[this_i] + min_clock_cycle*1E9/160.316E6:
            continue
        sel_hits.append(hit)

    for this_h in h:
        del this_h
    if not return_tref:
        return sel_hits
    else:
        return sel_hits, peakTime


import bisect
def getSumDensity(hit_collection, return_position=False, return_2ndhighest=False, return_hv=False, xlims=None, ylims=None, digi=True):

    x = [[] for i in range(10)]

    if not digi:
        fibres_seen = []
    
    
    for hit_i in hit_collection:
        if digi:
            this_station = hit_i.GetStation()-1
            this_orientation = hit_i.isVertical()

            scifiDet.GetSiPMPosition(hit_i.GetDetectorID(), a, b)
        else:
            if hit_i.GetDetectorID() in fibres_seen:
                continue
            else:
                fibres_seen.append(hit_i.GetDetectorID())
                
            this_station = hit_i.station()-1
            this_orientation = hit_i.orientation()

            if hit_i.row() != 1:
                continue
            
            scifiDet.GetPosition(hit_i.GetDetectorID(), a, b)
            
        
        if (xlims is not None) and this_orientation :
            if (a.X() + b.X())/2. < xlims[0]:
                continue
            if (a.X() + b.X())/2. >= xlims[1]:
                continue

        if (ylims is not None) and (not this_orientation) :
            if (a.Y() + b.Y())/2. < ylims[0]:
                continue
            if (a.Y() + b.Y())/2. >= ylims[1]:
                continue
        
        if this_orientation:
            x_i = (a.X()+b.X())/2.
        else:
            x_i = (a.Y()+b.Y())/2.

        bisect.insort(x[this_orientation*5+this_station], x_i)

    if return_position:
        pos_density_weighted = [None for i in range(10)]
        
    sum_density = 0
    sum_density_ver = 0
    sum_density_hor = 0
    if return_2ndhighest:
        sum_density2 = 0
        sum_density_ver2 = 0
        sum_density_hor2 = 0
    for i_plane in range(5):
        this_sum_density = 0
        this_sum_density_ver = 0
        this_sum_density_hor = 0
        for orientation in [0, 1]:
            if return_position:
                this_pos_density_weighted = 0
                this_sum_density_orientation = 0
            this_i = i_plane + orientation*5
            for i_hit, hit in enumerate(x[this_i]):
                
                this_weight = 0
                # Go backwards:
                start_range = i_hit-1
                if return_position:
                    start_range = i_hit
                for j_hit in range(start_range, -1, -1):
                    if (hit - x[this_i][j_hit]) < 1.:
                        this_weight += 1
                    else:
                        break
                # Now go forwards:
                for j_hit in range(i_hit+1, len(x[this_i])):
                    if (x[this_i][j_hit] - hit) < 1.:
                        this_weight += 1
                    else:
                        break
                if return_position:
                    this_pos_density_weighted += hit*this_weight
                    this_sum_density_orientation += this_weight
                this_sum_density += this_weight
                if orientation:
                    this_sum_density_ver += this_weight
                else:
                    this_sum_density_hor += this_weight
            if return_position:
                if this_sum_density_orientation:
                    pos_density_weighted[this_i] = this_pos_density_weighted/this_sum_density_orientation
        if this_sum_density > sum_density:
            if return_2ndhighest:
                sum_density2 = sum_density
                sum_density_ver2 = sum_density_ver
                sum_density_hor2 = sum_density_hor
            sum_density = this_sum_density
            sum_density_ver = this_sum_density_ver
            sum_density_hor = this_sum_density_hor
        if return_2ndhighest:
            if this_sum_density < sum_density and this_sum_density > sum_density2:
                sum_density2 = this_sum_density
                sum_density_ver2 = this_sum_density_ver
                sum_density_hor2 = this_sum_density_hor

                
    if return_position:
        return pos_density_weighted

    ret = []
    ret.append(min(sum_density_ver, sum_density_hor))
    
    if return_2ndhighest:
        ret.append(min(sum_density_ver2, sum_density_hor2)) # Return minimum density of plane with highest average density
        if return_hv:
            ret.append(sum_density_ver)
            ret.append(sum_density_hor) # Return minimum density of plane with highest average density

    return (*ret,)

def getSciFiAngle(hits):
    
    gv = ROOT.TGraph()
    gh = ROOT.TGraph()

    for hit in hits:
        if not hit.isValid():
            continue

        scifiDet.GetSiPMPosition(hit.GetDetectorID(), a, b)

        if hit.isVertical():
            gv.SetPoint(gv.GetN(), (a.Z()+b.Z())/2., (a.X()+b.X())/2.)
        else:
            gh.SetPoint(gh.GetN(), (a.Z()+b.Z())/2., (a.Y()+b.Y())/2.)
                
    fitv = gv.Fit("pol1", "SQN")
    fith = gh.Fit("pol1", "SQN")

    slopev = fitv.Get().GetParams()[1]
    interceptv = fitv.Get().GetParams()[0]

#    reducedchi2v = fitv.Get().Chi2()/fitv.Get().Ndf()
    reducedchi2v = fitv.Get().Chi2()

    slopeh = fith.Get().GetParams()[1]
    intercepth = fith.Get().GetParams()[0]
    
    reducedchi2h = fith.Get().Chi2()/fith.Get().Ndf()
    reducedchi2h = fith.Get().Chi2()


    reducedchi2both = (fitv.Get().Chi2() + fith.Get().Chi2())/(fitv.Get().Ndf()+fith.Get().Ndf())
    
    return slopev, slopeh, reducedchi2v, reducedchi2h, reducedchi2both

def tdiff(scifi_hits, mufilter_hits, MC = False):

    if MC:
        TDC2ns = 1
    else:
        TDC2ns = 1E9/160.316E6
    
    min_scifi_t = 1e9
    for hit in scifi_hits:
        if not hit.isValid():
            continue
        
        time = hit.GetTime()*TDC2ns
        if not MC:
            time = scifiDet.GetCorrectedTime(hit.GetDetectorID(), time,0)
        if time < min_scifi_t:
            min_scifi_t = time

    max_US_t = -1e9
    for hit in mufilter_hits:
        if not hit.isValid():
            continue
        if not hit.GetSystem() == 2:
            continue
        avg_hit_time = 0
        n_hits = 0
        for i in range(16):
            if i in [2, 5, 10, 13]:
                # Skip small SiPMs
                continue
            if hit.GetSignal(i) > -998:
                avg_hit_time += hit.GetTime(i)
                n_hits += 1
        if avg_hit_time/n_hits > max_US_t:
            max_US_t = avg_hit_time/n_hits

    return max_US_t - min_scifi_t

#def getSciFiAngle(hits):
#    
#    z_planes = [298.97, 311.97, 324.97, 337.97, 350.97]
#        
#    # calculate density-weighted average
#    positions = getSumDensity(hits, True)
#
#    gv = ROOT.TGraph()
#    gh = ROOT.TGraph()
#    
#    for i_station in range(5):
#        if positions[i_station]:
#            gv.SetPoint(gv.GetN(), z_planes[i_station], positions[i_station])
#        if positions[i_station+5]:
#            gh.SetPoint(gh.GetN(), z_planes[i_station], positions[i_station+5])
#            
#    fitv = gv.Fit("pol1", "SQN")
#    fith = gh.Fit("pol1", "SQN")
#
#    slopev = fitv.Get().GetParams()[1]
#    interceptv = fitv.Get().GetParams()[0]
#
#    reducedchi2v = fitv.Get().Chi2()
#
#    slopeh = fith.Get().GetParams()[1]
#    intercepth = fith.Get().GetParams()[0]
#    reducedchi2h = fith.Get().Chi2()
#
#
#    reducedchi2both = (fitv.Get().Chi2() + fith.Get().Chi2())
#    
#    return slopev, slopeh, reducedchi2v, reducedchi2h, reducedchi2both
