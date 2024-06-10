#!/usr/bin/env python3
from shapely.geometry import LineString

import simcem

db=simcem.defaultDatabase()

class InfeasibleDesign(Exception):
    pass

class Line():
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        self.m = (p2[1] - p1[1]) / (p2[0] - p1[0])
        self.c = p1[1] - self.m * p1[0]

    def intersection(l):
        xdiff = (self.p1[0] - self.p2[0], l.p1[0] - l.p2[0])
        ydiff = (self.p1[1] - self.p2[1], l.p1[1] - l.p2[1])

        def det(a, b):
            return a[0] * b[1] - a[1] * b[0]

        div = det(xdiff, ydiff)
        if div == 0:
            return None, None

        d = (det(self.p1, self.p2), det(l.p1, l.p2))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        return x, y

    def __call__(self, x):
        return self.m * x + self.c

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        return None,None

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

def full_qline(q, xf):
    # Calculate the q line first
    if abs(q-1) < 1e-5:
        #For q close to one, use differences in y to generate the two line
        #points (to avoid divide by zero).
        y1 = 1.0
        x1 = (y1 * (q-1) + xf) / q
        y2 = 0.0
        x2 = (y2 * (q-1) + xf) / q
    else:
        #For q away from 1, use differences in x to generate the two line
        #points (to avoid divide by zero).
        x1 = 1.0
        y1 = (x1 * q - xf) / (q-1)
        x2 = 0.0
        y2 = (x2 * q - xf) / (q-1)
    return LineString([(x1,y1), (x2,y2)])

def buildLines(xw, xf, xd, q, R):
    #Calculate the full q line
    q_line = full_qline(q, xf)
    #Calculate the full upper operating line
    upper_op_line = LineString([(0,xd / (R + 1)), (xd, xd)])
    #Find what is the crossover between the q line and operating line
    crossover = upper_op_line.intersection(q_line)
    if not crossover:
        raise InfeasibleDesign()
    #Now create the lower operating line
    lower_op_line = LineString([(xw,xw), crossover.coords[0]])
    #Trim the q and upper lines to the crossover intersection
    q_line = LineString([(xf, xf), crossover.coords[0]])
    upper_op_line = LineString([(xd, xd), crossover.coords[0]])
    
    return upper_op_line, lower_op_line, q_line, crossover

def minRfinder(xw, xf, xd, q, VLE):
    if xw > xf:
        raise RuntimeError("xw > xf?")
    if xf > xd:
        raise RuntimeError("xf > xd?")
    
    #First check infinte reflux will work
    inf_reflux_line = LineString([(xw,xw), (xd,xd)])
    if inf_reflux_line.intersection(VLE):
        #There's an intersection of the VLE line with the 45 degree
        #line, possible azeotrope in the way therefore there is no
        #minimum reflux that can get the job done! Not even infinite
        #will work (so we can't return that)
        raise InfeasibleDesign()

    #Next, just test if a flash drum would get it done (q-line/VLE intersection)
    full_q_line = full_qline(q, xf)
    i = full_q_line.intersection(VLE)

    if not i:
        raise RuntimeError("q line missed the VLE? Check your VLE data doesn't have gaps?")

    if (i.coords[0][0] < xw) and (i.coords[0][1] > xd):
        #Yup, a flash drum would do it. No reflux needed.
        return 0
    
    #We can still have a solution for R=0 if its a stripping column
    #(feed is top tray and provides the liquid). Check if 0 reflux is
    #a valid solution.
    try:
        upper_op_line, lower_op_line, q_line, crossover = buildLines(xw=xw, xf=xf, xd=xd, q=q, R=0)
        if not upper_op_line.intersection(VLE) and not lower_op_line.intersection(VLE):
            #The R=0 design does not intersect the VLE, so its fine to use!
            return 0
    except InfeasibleDesign as e:
        #If the R=0 design is not feasible, just carry on with the bracketing below
        pass
    
    #Guess a reflux ratio, then exponentially try to bracket the good/bad limits
    R1 = 1.0
    try:
        upper_op_line, lower_op_line, q_line, crossover = buildLines(xw=xw, xf=xf, xd=xd, q=q, R=R1)
        #Test if the op lines intersect the VLE lines
        R1_bad = bool(upper_op_line.intersection(VLE)) or bool(lower_op_line.intersection(VLE))
    except InfeasibleDesign as e:
        #If we can't even build the op/q lines in the graph then its definitely "bad"
        R1_bad = True
        
    R2 = R1
    while True:
        if R1_bad:
            R2 *= 2
        else:
            R2 /= 2

        try:
            upper_op_line, lower_op_line, q_line, crossover = buildLines(xw=xw, xf=xf, xd=xd, q=q, R=R2)
            R2_bad = bool(upper_op_line.intersection(VLE)) or bool(lower_op_line.intersection(VLE))
        except InfeasibleDesign as e:
            #If we can't even build the op/q lines in the graph then its definitely "bad"
            R2_bad = True

        if R2_bad != R1_bad:
            break
        R1 = R2
        
    if R1_bad:
        Rlow = R1
        Rhigh = R2
    else:
        Rlow = R2
        Rhigh = R1

    while (Rhigh - Rlow)/ Rhigh > 0.0001:
        Rmid = (Rhigh + Rlow) / 2
        try:
            upper_op_line, lower_op_line, q_line, crossover = buildLines(xw=xw, xf=xf, xd=xd, q=q, R=Rmid)
            if upper_op_line.intersection(VLE) or lower_op_line.intersection(VLE):
                Rlow = Rmid
            else:
                Rhigh = Rmid
        except InfeasibleDesign as e:
            #If we can't even build the op/q lines in the graph then its definitely "bad"
            Rlow = Rmid

    return Rhigh # Rlow is by definition, bad. So return something that is known good

def RPlot(ax, xd, R, xw=None, xf=None, q=None):
    '''Add a line for the correspondance of the upper operating line and the reflux ratio'''
    ax.plot([0, xd], [xd / (R+1), xd], color="teal")
    ax.text(0, xd / (R+1), f"{xd / (R+1):.3f}", color="teal", horizontalalignment="right", verticalalignment="center")
    if xf and xw and q:
        upper_op_line, lower_op_line, q_line, crossover = buildLines(xw, xf, xd, q, R)
        ax.plot([x for x,y in lower_op_line.coords], [y for x,y in lower_op_line.coords], 'teal', )

def minStages(xw, xd, VLE):
    inf_reflux_line = LineString([(xw,xw), (xd,xd)])

    stages = stepping(inf_reflux_line, VLE, xstart=xw)
    
    if len(stages) > 1:
        final_frac = (xd - stages[-2][1]) / (stages[-1][1] - stages[-2][1])
    else:
        final_frac = (xd - xw) / (stages[-1][1] - xw)

    return stages, len(stages)-1 + final_frac


def stepping(op_line, VLE, xstart=None, ystart=None, up=True, max_stages=300, murphree=lambda x:1):
    if xstart is None and ystart is None:
        raise RuntimeError("Must start the stepping somewhere?")
        
    #print('####', op_line)
    stages=[]
    if up:
        #We're stepping upwards
        if xstart is None:
            #We've been given a y value, start the stepping with the first intersection
            i = LineString([(0,ystart), (1,ystart)]).intersection(op_line)
            if not i:
                return []
            #All seemed to work, so carry on
            xstage = i.coords[0][0]
            yop = i.coords[0][1]
        else:
            xstage = xstart
            i = LineString([(xstart,0), (xstart,1)]).intersection(op_line)
            yop = i.coords[0][1]
            
        #Limit the number of stages incase it gets stuck
        for idx in range(max_stages):
            #Calculate the stage's y value at the VLE intersection
            ystage_equil = LineString([(xstage,0), (xstage,1)]).intersection(VLE).coords[0][1]
            eff = murphree(xstage)
            ystage = ystage_equil * eff + yop * (1 - eff)
            stages.append((xstage, ystage))
            yop=ystage
            #Calculate the next intersection on the operating line
            i = LineString([(0,ystage), (1,ystage)]).intersection(op_line)
            #If we're past the end of the operating line, stop. Same check
            #as before for no intersection.
            if not i:
                return stages
            
            # We've not past it, so carry on
            xstage = i.coords[0][0]
    else:
        #We're stepping downwards
        if ystart is None:
            #We've been given a y value, start the stepping with the first intersection
            i = LineString([(xstart,0), (xstart,1)]).intersection(op_line)
            if not i:
                return []
            #All seemed to work, so carry on
            xstage = i.coords[0][0]
            yop = i.coords[0][1]
        else:
            xstage = xstart
            i = LineString([(xstart,0), (xstart,1)]).intersection(op_line)
            yop = i.coords[0][1]
        for idx in range(max_stages):
            #Calculate the stage's x value at the VLE intersection
            i = LineString([(0,yop), (1,yop)]).intersection(VLE)
            ystage = i.coords[0][1]
            xstage = i.coords[0][0]
            stages.append((xstage, ystage))
            
            #Calculate the next intersection on the operating line
            i = LineString([(xstage,0), (xstage,1)]).intersection(op_line)

            #If we're past the end of the operating line, stop. Same check
            #as before for no intersection.
            if not i:
                return stages            
            yop = i.coords[0][1]

    #Reached the max number of stages
    raise InfeasibleDesign()

def mccabe_thiele(xw, xf, xd, q, R, VLE, murphree = lambda x,upper : 1, skipsolve=False, direction=None):
    if direction is None:
        direction='up'
    if xw > xf:
        raise RuntimeError("xw > xf?")
    if xf > xd:
        raise RuntimeError("xf > xd?")
    if R < 0:
        raise RuntimeError("R < 0?")

    #No need to catch InfeasibleDesign exceptions. Happy to let it propogate up the stack
    upper_op_line, lower_op_line, q_line, crossover = buildLines(xw, xf, xd, q, R)

    if skipsolve:
        return [], 0, 0, upper_op_line, lower_op_line, q_line

    #Now perform the stepping for the upper operating line
    if direction == 'up':
        low_stages = stepping(lower_op_line, VLE, xstart=xw, murphree=lambda x: murphree(x, False), up=True)
        upper_stages = stepping(upper_op_line, VLE, ystart=low_stages[-1][1], murphree=lambda x: murphree(x, True), up=True)
    elif direction == "down":
        upper_stages = stepping(upper_op_line, VLE, xstart=xd, murphree=lambda x: murphree(x, True), up=False)
        low_stages = stepping(lower_op_line, VLE, xstart=upper_stages[-1][0], murphree=lambda x: murphree(x, False), up=False)
        low_stages.reverse()
        upper_stages.reverse()
        #Our downward stepping is too greedy for the "feed" tray calcs
        low_stages.append(upper_stages[0])
        upper_stages = upper_stages[1:]
    elif direction == "mid":
        start = full_qline(q, xf).intersection(VLE)
        upper_stages = stepping(upper_op_line, VLE, ystart=start.coords[0][1], murphree=lambda x: murphree(x, True), up=True)
        low_stages = stepping(lower_op_line, VLE, xstart=start.coords[0][0], murphree=lambda x: murphree(x, False), up=False)
        low_stages.reverse()
        low_stages.append((start.coords[0][0], start.coords[0][1]))     
    else:
        raise RuntimeError('Unrecognised stepping direction')    
    
    feed_stage_ID = len(low_stages)
    stages = low_stages + upper_stages

    #The stages[-2][1] != stages[-1][1] check is to prevent situations
    #where the stepping has gotten stuck from causing a div by zero
    if len(stages) > 1 and stages[-2][1] != stages[-1][1]:
        final_frac = (xd - stages[-2][1]) / (stages[-1][1] - stages[-2][1]) + (stages[1][0] - xw) / (stages[1][0] - stages[0][0])
    else:
        final_frac = (xd - xw) / (stages[-1][1] - xw)

    stage_count = len(stages) - 2 + final_frac

    return stages, stage_count, feed_stage_ID, upper_op_line, lower_op_line, q_line

def xyPlot(VLE, xf=None, xw=None, xd=None, stages=None, q_line=None, upper_op_line=None, lower_op_line=None, yx=True):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
    fig,ax = plt.subplots(1, 1)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_xlim((0, 1)); ax.set_ylim((0, 1))
    ax.set_aspect('equal')
    ax.minorticks_on()
    fig.gca().xaxis.set_minor_locator(MultipleLocator(0.02))
    fig.gca().yaxis.set_minor_locator(MultipleLocator(0.02))
    fig.gca().xaxis.set_major_locator(MultipleLocator(0.1))
    fig.gca().yaxis.set_major_locator(MultipleLocator(0.1))
    ax.grid(which='minor', axis='both', color='silver', linestyle='-')
    ax.grid(which='major', axis='both', color='dimgray', linestyle='-')
    if yx:
        ax.plot([0,1], [0,1], color='black')  # y=x
    ax.plot([x for x,y in VLE.coords], [y for x,y in VLE.coords], 'b')  # curve

    if xf:
        ax.plot([xf,xf],[0,xf], color='green', linestyle='dashed')  # xw
        ax.text(xf, 0, '$x_F$', color='green')
    if xw:
        ax.plot([xw,xw],[0,xw], color='green', linestyle='dashed')  # xw
        ax.text(xw, 0, '$x_W$', color='green')
    if xd:
        ax.plot([xd,xd],[0,xd], color='green', linestyle='dashed')  # xw
        ax.text(xd, 0, '$x_D$', color='green')
    if q_line:
        ax.plot([x for x,y in q_line.coords], [y for x,y in q_line.coords], color='red')
    if upper_op_line:
        ax.plot([x for x,y in upper_op_line.coords], [y for x,y in upper_op_line.coords], color='teal')
    if lower_op_line:
        ax.plot([x for x,y in lower_op_line.coords], [y for x,y in lower_op_line.coords], color='teal')
    if stages:
        yold = xw
        stepping_points = []
        for xstage, ystage in stages:
            stepping_points.append((xstage,yold))
            stepping_points.append((xstage,ystage))
            yold = ystage
        ax.plot([x for x,y in stepping_points], [y for x,y in stepping_points], 'k-')
    return fig, ax

import math
import numpy

def solveVLE(activityA, activityB, vapourPA, vapourPB, targetP=760, minx=0, maxx=1.0, stepx=0.01,
             Thigh_start=120+273.15, Tlow_start=30+273.15):
    def PartialP(x, T, activity, vapourP):
        return x * vapourP(T) * activity(x)

    T=[]
    x=[]
    y=[]
    for xval in numpy.arange(minx, maxx, stepx):
        Thigh=Thigh_start
        Tlow=Tlow_start
        while ((Thigh-Tlow) > 0.001):
            Tmid=0.5*(Thigh+Tlow)
            PA=PartialP(xval, Tmid, activityA, vapourPA)
            PB=PartialP(1-xval, Tmid, activityB, vapourPB)
            if ((PA+PB)>targetP):
                Thigh=Tmid
            else:
                Tlow=Tmid
        Tval = 0.5*(Thigh+Tlow)
        T.append(Tval)
        x.append(xval)
        y.append(PartialP(xval, Tval, activityA, vapourPA)/targetP)
    return T,x,y

def ActivityMargules(x, AB, BA):
    return math.exp((1-x)*(1-x)*(AB + 2 * x * (BA-AB)))

def ActivityvanLaar(x, AB, BA):
    return math.exp(AB*((BA*(1-x)/(AB*x+BA*(1-x)))**2))

def IPAVLE(bar=1, mode=1):
    #The torr data and solve commands are from the old resource, but this function now uses much better data from https://doi.org/10.1021/je9503113

    def watertorr(T):
        T = T-273.15
        return 10**(8.07131-1730.630/(T+233.426))
    
    def IPAtorr(T):
        T = T-273.15
        return 10**(8.87829-2010.320/(T+252.636))
    
    def torrToBar(P):
        return (P/760)*1.01325
        
    def ActivityvanLaarIPA(x):
        return ActivityvanLaar(x, 2.4702, 1.0938)
    
    def ActivityvanLaarWater(x):
        return ActivityvanLaar(x, 1.0938, 2.4702)
    
    def ActivityMargulesIPA(x):
        return ActivityMargules(x, 2.3319, 0.8976)
    
    def ActivityMargulesWater(x):
        return ActivityMargules(x, 0.8976, 2.3319)
    
    def ActivityIdeal(x):
        return 1
    
    def psatWater(T):
        A=+16.5700
        B=+3984.92
        C=-39.724
        return math.exp(A - B / (C+T))*1e3

    def psatIPA(T):
        A=+16.4089
        B=+3439.60
        C=-63.417
        return math.exp(A - B / (C+T))*1e3

    #Mode 0 and 1 are the paper, mode 2 and 3 are Perrys? and don't catch the azeotrope right in T
    if mode == 0:
        return solveVLE(lambda x : ActivityMargules(x, 2.2101, 1.0443), lambda x : ActivityMargules(x, 1.0443, 2.2101), psatIPA, psatWater, 1e5 * bar, 0.0, 1.0011, 0.001, Tlow_start = 200, Thigh_start = 500)
    elif mode == 1:
        return solveVLE(lambda x : ActivityvanLaar(x, 2.3108, 1.1810), lambda x : ActivityvanLaar(x, 1.1810, 2.3108), psatIPA, psatWater, 1e5 * bar, 0.0, 1.0011, 0.001, Tlow_start = 200, Thigh_start = 500)
    elif mode == 2:
        return solveVLE(ActivityMargulesIPA, ActivityMargulesWater, IPAtorr, watertorr, 760 * bar, 0.0, 1.0011, 0.001, Tlow_start = 200, Thigh_start = 500)
    elif mode == 3:
        return solveVLE(ActivityvanLaarIPA, ActivityvanLaarWater, IPAtorr, watertorr, 760 * bar, 0.0, 1.0011, 0.001, Tlow_start = 200, Thigh_start = 500)
    else:
        raise RuntimeError("Bad mode")

def benzeneVapP(T):
    return 10**(4.72583 - 1660.652 / (T - 1.461))

#https://webbook.nist.gov/cgi/cbook.cgi?ID=C108883&Units=SI&Mask=4#Thermo-Phase
def tolueneVapP(T):
    return 10**(4.14157 - 1377.578 / (T - 50.507))

def heptaneVapP(T):
    return 10**(4.02832 - 1268.636 / (T - 56.199))

def benzeneTolueneVLE(bar=1):
    return solveVLE(lambda x: 1, lambda x: 1, benzeneVapP, tolueneVapP, bar, 0.0, 1.0011, 0.001, Tlow_start = 200, Thigh_start = 500)

def benzeneHeptaneVLE(bar=1):
    return solveVLE(lambda x: 1, lambda x: 1, benzeneVapP, heptaneVapP, bar, 0.0, 1.0011, 0.001, Tlow_start = 200, Thigh_start = 500)

def ethanolWaterVLE(bar=1):
    #https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=4#Thermo-Phase
    def waterVapP(T):
        return 10**(4.6543 - 1435.264 / (T - 64.848))
    #https://webbook.nist.gov/cgi/cbook.cgi?Name=ethanol&Units=SI
    def ethVapP(T):
        return 10**(5.37229 - 1670.409 / (T - 40.191))
    
    def ActivityMargulesEth(x):
        return ActivityMargules(x, 1.6022, 0.7947)
    def ActivityMargulesWater(x):
        return ActivityMargules(x, 0.7947, 1.6022)
    def ActivityvanLaarEth(x):
        return ActivityvanLaar(x, 1.6798, 0.9227)
    def ActivityvanLaarWater(x):
        return ActivityvanLaar(x, 0.9227, 1.6798)

    return solveVLE(ActivityMargulesEth, ActivityMargulesWater, ethVapP, waterVapP, bar, 0, 1.011, 0.001, Tlow_start=350, Thigh_start=500)

def ethanolWaterH(Ts, xs, bar=1):
    #Mixing enthalpies taken from here
    #http://dx.doi.org/10.1016/0021-9614(75)90261-X
    def coeffs(T):
        C={0.0:[-3.63868, +1.83829, -2.32763], 0.5: [+9.25982, -4.83586, +6.37228],
           1.5:[-14.04894, +7.51661, -10.11280], 2.5: [10.91318, -5.89498, +7.98868],
           4.5:[-2.79986, +1.50557, -2.03127]}
        
        return [(n, bcd[0] * 1e5 + bcd[1] * 1e3 * T + bcd[2] * T * T) for n,bcd in C.items()]
    
    B=[1580.0, 1785.0, 3487.0, 3187.0, 1957.0]
    def mixEnthalpy(x, T):
        return 0
        C=coeffs(T)
        return x * (1.0 - x) * sum([a * (x**n) for n,a in C])

    hwaterform = simcem.ModelIncompressibleLiquid(db, simcem.Components({'H2O':1}), p=1e5*bar, T=298.15).H()
    hethform = simcem.ModelIncompressibleLiquid(db, simcem.Components({'C2H5OH':1}), p=1e5*bar, T=298.15).H()
    href = simcem.ModelIncompressibleLiquid(db, simcem.Components({'H2O':1}), p=1e5*bar, T=373.15).H() - hwaterform
    Hliq = []
    Hvap = []
    for T,x in zip(Ts, xs):
        liq = simcem.ModelIncompressibleLiquid(db, simcem.Components({'H2O':1-x, 'C2H5OH':x}), p=1e5*bar, T=T)
        gas = simcem.ModelIdealGasTp(db, simcem.Components({'H2O':1-x, 'C2H5OH':x}), p=1e5*bar, T=T)
        Hliq.append(liq.H() - x * hethform - (1.0 - x) * hwaterform - href + mixEnthalpy(x, T))
        Hvap.append(gas.H() - x * hethform - (1.0 - x) * hwaterform - href)
    return Hliq, Hvap

def vapourEnthalpy(x, T):
    heth = speciesData["C2H5OH"].Hf0(T, 'Gas') - href - hethform
    hwater = speciesData["H2O"].Hf0(T, 'Gas') - href - hwaterform
    return x * heth + (1.0 - x) * hwater

def pentaneHexaneVLE():
    #1 atm, taken from this example
    #https://www.cpp.edu/~tknguyen/che313/pdf/chap4-4c.pdf
    data = [342.06, 0.00000, 0.00000,
            339.40, 0.05000, 0.12705,
            336.91, 0.10000, 0.23699,
            334.58, 0.15000, 0.33263,
            332.39, 0.20000, 0.41626,
            330.32, 0.25000, 0.48975,
            328.38, 0.30000, 0.55462,
            326.53, 0.35000, 0.61214,
            324.79, 0.40000, 0.66335,
            323.14, 0.45000, 0.70911,
            321.56, 0.50000, 0.75016,
            320.07, 0.55000, 0.78711,
            318.64, 0.60000, 0.82048,
            317.28, 0.65000, 0.85070,
            315.97, 0.70000, 0.87816,
            314.72, 0.75000, 0.90317,
            313.53, 0.80000, 0.92601,
            312.38, 0.85000, 0.94692,
            311.28, 0.90000, 0.96610,
            310.22, 0.95000, 0.98374,
            309.20, 1.00000, 1.00000]
    return data[0::3], data[1::3], data[2::3]


def HxyPlot(x=None, y=None, Hliq=None, Hvap=None, xf=None, hf=None, xd=None, xw=None, tielinecount=None, stages=[], oplines=[], hdmin=None, hwmin=None, hwprime=None, hdprime=None):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
    fig,ax = plt.subplots(1, 1)
    ax.set_xlabel('$x$')
    ax.set_ylabel(r'$h~\left(\mathrm{kJ}~\mathrm{kmol}^{-1}\right)$')
    ax.set_xlim((0, 1)); #ax.set_ylim((0, 1)); ax.set_aspect('equal')
    ax.minorticks_on()
    fig.gca().xaxis.set_minor_locator(MultipleLocator(0.02))
    fig.gca().xaxis.set_major_locator(MultipleLocator(0.1))
    ax.grid(which='minor', axis='both', color='silver', linestyle='-')
    ax.grid(which='major', axis='both', color='dimgray', linestyle='-')
    if x is not None and Hliq is not None:
        ax.plot(x, Hliq, 'b')
    if x is not None and Hvap is not None:
        ax.plot(x, Hvap, 'r')

    if x is not None and y is not None and Hliq is not None and Hvap is not None and tielinecount is not None:
        import scipy.interpolate
        fhv = scipy.interpolate.interp1d(x, Hvap, kind='cubic')
        step = int(len(x) / tielinecount)
        for xv, yv, hl in zip(x[::step], y[::step], Hliq[::step]):
            ax.plot([xv, yv], [hl, fhv(yv)], 'b-')
        
    for stage in oplines:
        ax.plot([stage[0][0], stage[1][0]], [stage[0][1], stage[1][1]], '-', color="teal")

    for stage in stages:
        ax.plot([stage[0][0], stage[1][0]], [stage[0][1], stage[1][1]], 'b-x')

    if hdmin is not None and xd is not None:
        ax.plot([xd], [hdmin], '.', color='dimgray')
        ax.text(xd, hdmin, '$P_C^{(min)}$', fontsize=12, color='dimgrey')

    if hwmin is not None and xw is not None:
        ax.plot([xw], [hwmin], '.', color='dimgray')
        ax.text(xw, hwmin, '$P_R^{(min)}$', fontsize=12, color='dimgrey')

    if xw is not None and xf is not None and xd is not None and hwprime is not None and hf is not None and hdprime is not None:
        ax.plot((xw, xf, xd), (hwprime, hf, hdprime), '.g-')

    # We have to save the current scaling, as it seems to be lost on autoscale???
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.autoscale(False, axis='both')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    
    if x is not None and Hliq is not None and y is not None and Hvap is not None:
        import scipy.interpolate
        fhl = scipy.interpolate.interp1d(x, Hliq, kind='cubic')
        fhv = scipy.interpolate.interp1d(x, Hvap, kind='cubic')
        bot, top = ax.get_ylim()
        for x, txt in [(xf, 'x_F'), (xd, 'x_D') , (xw, 'x_W')]:
            if x is not None:
                ax.plot([x, x], [bot, top], '--', color='black')
                ax.text(x, bot, '$'+txt+'$', fontsize=12, color='black')
        
    return fig, ax

def minRPS(x, y,  Hliq, Hvap, xf, hf, xw, xd, tielinecount=None):
    import scipy.interpolate
    fhv = scipy.interpolate.interp1d(x, Hvap, kind='cubic')
    fhl = scipy.interpolate.interp1d(x, Hliq, kind='cubic')
    hdprime = -1e30
    hwprime = +1e30
    if tielinecount is None:
        tielinecount = len(x)
    step = int(len(x) / tielinecount)
    ties = []
    for xv,yv,hl, in zip(x[::step],y[::step],Hliq[::step]):
        if (xv < xw and yv < xw) or (xv > xd and yv > xd):
            #Skip tielines completely past the bottoms/top lines
            continue
        hv = fhv(yv)
        tieline = ((xv,hl),(yv,hv))
        feedline = ((xf, -1), (xf, 1))
        fx,fh = line_intersection(tieline, feedline)
        if fx is None:
            raise InfeasibleDesign()

        if fh > hf:
            wx,wh = line_intersection(tieline, ((xw, -1), (xw,1)))
            ties.append(((yv,hv), (wx,wh)))
            hwprime = min(hwprime, wh)
        else:
            dx,dh = line_intersection(tieline, ((xd, -1), (xd,1)))
            ties.append(((xv,hl), (dx,dh)))
            hdprime = max(hdprime, dh)
    
    #Now check where the bottoms intersection is on the xd line via the feed point
    wx,wh = line_intersection(((xw, hwprime), (xf, hf)), ((xd, -1), (xd, 1)))
    hdmin = max(hdprime, wh)
    #ax.plot([xw], [hwprime], 'xg')
    #ax.plot([xd], [hdprime], 'xg')

    #The max is required here as graphically negative R values are
    #possible, but not physical
    minR = (hdmin - fhv(xd)) / (fhv(xd) - fhl(xd))
    if minR < 0:
        minR = 0
    return minR, ties, hwprime, hdprime, hdmin

    
def PonchonSavarit_down(x, y, Hliq, Hvap, xw, xd, xf, hf, R):
    import scipy.interpolate
    fhv = scipy.interpolate.interp1d(x, Hvap, kind='cubic')
    fhl = scipy.interpolate.interp1d(x, Hliq, kind='cubic')
    fVLE = scipy.interpolate.interp1d(x, y, kind='cubic')
    fELV = scipy.interpolate.interp1d(y, x, kind='cubic')
    
    hdprime = R * (fhv(xd) - fhl(xd)) + fhv(xd)
    mainline = ((xd, hdprime), (xf, hf))
    _, hwprime = line_intersection(mainline, ((xw, -1), (xw, 1)))
    mainline = ((xd, hdprime), (xw, hwprime))
    mainlinef = scipy.interpolate.interp1d((xd, xw), (hdprime, hwprime), kind='linear')
    hvline = LineString([(a, b) for a, b in zip(x, Hvap)])
    hlline = LineString([(a, b) for a, b in zip(x, Hliq)])

    oplines = []
    upper_stages = []
    xstage = xd
    hlstage = fhl(xstage)
    while True:
        #Create operating line to find y
        op_line = LineString([(xstage,hlstage), (xd,hdprime)])
        oplines.append([(xstage,hlstage), (xd,hdprime)])
        ystage, hvstage = op_line.intersection(hvline).coords[0]
        #Now use VLE data to find x again
        #print(ystage, hvstage)
        xstage = fELV(ystage)
        hlstage = fhl(xstage)
        if hlstage > mainlinef(xstage):
            break
        upper_stages.append([(xstage, hlstage), (ystage, hvstage)])

    lower_stages = [[(xstage, hlstage), (ystage, hvstage)]]
    
    while True:
        #Create operating line to find y
        l = Line((xw,hwprime), (xstage,hlstage))
        op_line = LineString([(xw,hwprime), (1, l(1))])
        ystage, hvstage = op_line.intersection(hvline).coords[0]
        oplines.append([(xw,hwprime), (ystage, hvstage)])
        
        #Now use VLE data to find x again
        xstage = fELV(ystage)
        hlstage = fhl(xstage)
        lower_stages.append([(xstage, hlstage), (ystage, hvstage)])
        if xstage < xw:
            break

    frac_stage = (lower_stages[-1][1][0] - xw) / (lower_stages[-1][1][0] - lower_stages[-1][0][0])
    stage_count = len(lower_stages) + len(upper_stages) - 1 + frac_stage
    return hdprime, hwprime, upper_stages, lower_stages, oplines, stage_count

def PonchonSavarit_up(x, y, Hliq, Hvap, xw, xd, xf, hf, R):
    import scipy.interpolate
    fhv = scipy.interpolate.interp1d(x, Hvap, kind='cubic')
    fhl = scipy.interpolate.interp1d(x, Hliq, kind='cubic')
    fVLE = scipy.interpolate.interp1d(x, y, kind='cubic')
    
    hdprime = R * (fhv(xd) - fhl(xd)) + fhv(xd)
    mainline = ((xd, hdprime), (xf, hf))
    _, hwprime = line_intersection(mainline, ((xw, -1), (xw, 1)))
    mainline = ((xd, hdprime), (xw, hwprime))
    mainlinef = scipy.interpolate.interp1d((xd, xw), (hdprime, hwprime), kind='linear')
    hvline = LineString([(a, b) for a, b in zip(x, Hvap)])
    hlline = LineString([(a, b) for a, b in zip(x, Hliq)])

    oplines = []
    lower_stages = []
    ystage = xw
    hvstage = fhv(ystage)
    for i in range(100):
        #Create operating line to find y
        op_line = LineString([(xw, hwprime), (ystage, hvstage)])
        oplines.append([(xw,hwprime), (ystage, hvstage)])
        
        xstage, hlstage = op_line.intersection(hlline).coords[0]
        ystage = fVLE(xstage)
        hvstage = fhv(ystage)
        lower_stages.append([(xstage, hlstage), (ystage, hvstage)])

        #First if condition guards against the feed stage somehow
        #going beyond xd the mainline range in one step
        if ystage > xd or hvstage <= mainlinef(ystage):
            break
        
    upper_stages=[]
    
    if ystage <= xd and hvstage > mainlinef(ystage):
        raise InfeasibleDesign()

    for i in range(100):
        if ystage >= xd:
            break
        #Create an operating line that we can extend to the Hliq line
        l = Line((xd,hdprime), (ystage, hvstage))
        op_line = LineString([(xd, hdprime), (0, l(0))])
        xstage, hlstage = op_line.intersection(hlline).coords[0]
        oplines.append([(xd,hdprime), (xstage, hlstage)])
        
        #Now use VLE data to find x again
        ystage = fVLE(xstage)
        hvstage = fhv(ystage)
        upper_stages.append([(xstage, hlstage), (ystage, hvstage)])

    if ystage < xd:
        raise InfeasibleDesign()

    stages = lower_stages + upper_stages
    frac_stage = (xd - stages[-1][0][0]) / (stages[-1][1][0] - stages[-1][0][0]) 
    stage_count = len(stages) - 1 + frac_stage
    return hdprime, hwprime, upper_stages, lower_stages, oplines, stage_count

import unittest
class TestDistillation(unittest.TestCase):

    
    def test_mccabe_thiele_design(self):
        #heptane-ethyl benzene design from EX3502 tutorials
        x = [0, 0.08, 0.250, 0.485, 0.790, 1.0]
        y = [0, 0.23, 0.514, 0.73, 0.904, 1.0]
        import scipy.interpolate 
        f = scipy.interpolate.interp1d(x, y, kind='cubic')
        import numpy as np
        x = [x for x in np.arange(0, 1.0, 0.01)]
        y = [f(x) for x in np.arange(0, 1.0, 0.01)]
        VLE = LineString(zip(x,y))
        xw = 0.011
        xf = 0.42
        xd = 0.97
        q = 0.5
        R = 2.5
        
        stages, stage_count, feed_stage_ID, upper_op_line, lower_op_line, q_line = mccabe_thiele(xw, xf, xd, q, R, VLE, skipsolve=False, direction='up')
        fig, ax = xyPlot(VLE, xf=xf, xw=xw, xd=xd, stages=stages, q_line=q_line, upper_op_line=upper_op_line, lower_op_line=lower_op_line)
        fig.savefig('up.png', bbox_inches='tight', pad_inches=0, dpi=500)
        self.assertAlmostEqual(stage_count, 11.88, delta=0.05)
        self.assertEqual(feed_stage_ID, 6)
        
        stages, stage_count, feed_stage_ID, upper_op_line, lower_op_line, q_line = mccabe_thiele(xw, xf, xd, q, R, VLE, skipsolve=False, direction='down')
        fig, ax = xyPlot(VLE, xf=xf, xw=xw, xd=xd, stages=stages, q_line=q_line, upper_op_line=upper_op_line, lower_op_line=lower_op_line)
        fig.savefig('down.png', bbox_inches='tight', pad_inches=0, dpi=500)
        self.assertAlmostEqual(stage_count, 11.88, delta=0.05)
        self.assertEqual(feed_stage_ID, 6)
        
        stages, stage_count, feed_stage_ID, upper_op_line, lower_op_line, q_line = mccabe_thiele(xw, xf, xd, q, R, VLE, skipsolve=False, direction='mid')
        fig, ax = xyPlot(VLE, xf=xf, xw=xw, xd=xd, stages=stages, q_line=q_line, upper_op_line=upper_op_line, lower_op_line=lower_op_line)
        fig.savefig('mid.png', bbox_inches='tight', pad_inches=0, dpi=500)
        self.assertAlmostEqual(stage_count, 11.88, delta=0.05)
        self.assertEqual(feed_stage_ID, 6)
        
    def test_mccabe_thiele_minR_minstages(self):
        #Class example from here
        #https://www.cpp.edu/~tknguyen/che313/pdf/chap4-4d.pdf
        T,x,y = pentaneHexaneVLE()
        VLE = LineString(zip(x,y))
        xw = 0.1
        xf = 0.40
        xd = 0.9
        q = 1
        Rmin = minRfinder(xw, xf, xd, q, VLE)
        minstages, minstagecount = minStages(xw, xd, VLE)
        self.assertAlmostEqual(Rmin, 0.9, delta=0.01)
        self.assertAlmostEqual(minstagecount, 4.1, delta=0.1)

    def test_ponchonsavarit(self):
        T,x,y = ethanolWaterVLE()
        Hliq, Hvap = ethanolWaterH(T, x, bar=1e5)
        xw = 0.1
        xd = 0.7
        xf = 0.2
        hf = 1000
        minR, ties, hwprime, hdprime, hdmin = minRPS(x, y,  Hliq, Hvap, xf, hf, xw, xd, tielinecount=20)
        R = 1.5 * minR
        self.assertAlmostEqual(minR, 0.4658471180244332, delta=0.01)
        
        hdprime, hwprime, upper_stages, lower_stages, oplines, stage_count = PonchonSavarit_up(x=x, y=y, Hliq=Hliq, Hvap=Hvap, xw=xw, xd=xd, xf=xf, hf=hf, R=R)
        self.assertAlmostEqual(stage_count, 6.563957204216382, delta=0.01)
        self.assertEqual(len(lower_stages), 2)
    
    # Future unit tests!
    #    #https://www.cpp.edu/~tknguyen/che313/pdf/chap4-4c.pdf
    #    #https://www.cpp.edu/~tknguyen/che313/pdf/chap4-5.pdf

if __name__ == '__main__':
    #T,x,y = benzeneTolueneVLE()
    #VLE = LineString(zip(x,y))
    #x= {'xw': 0.057, 'xf': 0.23, 'xd': 0.79, 'Rfactor': 1.18, 'q': 0.674,}# 'Rmin': 1.999969482421875}
    #x['Rmin'] = minRfinder(xw=x['xw'], xf=x['xf'], xd=x['xd'], q=x['q'], VLE=VLE)
    #x['R'] = x['Rmin'] * x['Rfactor']
    #stages, stage_count, feed_stage_ID, upper_op_line, lower_op_line, q_line = mccabe_thiele(x['xw'], x['xf'], x['xd'], x['q'], x['R'], VLE)
    #
    #fig, ax = xyPlot(VLE, xf=0.21, xw=0.059, xd=0.73, stages=stages, q_line=q_line, upper_op_line=upper_op_line, lower_op_line=lower_op_line)
    #fig.savefig('out.png', bbox_inches='tight', pad_inches=0, dpi=500)
    unittest.main()
