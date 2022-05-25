# IMPORTANT: This is _not_ an exercise in how to write nice Python
# INSTEAD this is an exercise in expressing a single algorithm in 3 languages
# Mathematica, R and Python
# Data structures are INTENTIONALLY _not_ 'Pythonic'
# Data structures (eg strings) have been chosen as 
# efficient and well supported _across languages_

# IF the same deterministic algorithm gives the same output on the same input
# _across multiple languages_
# THEN we are certain that language-specific behaviour has been eliminated
# (And that is a very powerful code property that is _much_ more important than
#  being fashionable in any particular language)

# For Python: numpy includes most of the array functionality the diem algorithm uses 

import numpy as np
from more_itertools import pairwise
from collections import Counter
import multiprocessing
import time
import os


#WARNING conversion to ERROR:
np.warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

# diem helper functions
# [direct replacements for Mathematica diem helper functions]
# function names are intended as self explanatory
#    clarifications in comments
# ________________________

def Total(listlikeArg):
    """Total([[1,2,3],[3,2,1]]) == array([4, 4, 4])"""
    return np.sum(np.array(listlikeArg),axis = 0)


def l2mers(listlikeArg):
    """The kmers of a list, where k==2"""
    return np.array([[a,b] for a,b in pairwise(listlikeArg)])

# The Mca version of this is more general. Only NaN insertion is needed by diem however
def MultiInsertnan(intoList, positions):
    outList = list(intoList) # copy
    for x in positions:
        outList[x:x] = np.nan
    return outList

#def StringReplace(text, charsFrom, charsTo):
#    return text.translate(str.maketrans(charsFrom, charsTo))

StringReplace20_dict = str.maketrans('02', '20')
def StringReplace20(text):
    """will _!simultaneously!_ replace 2->0 and 0->2"""
    return text.translate(StringReplace20_dict)

# EPDPfix:
def listOfNumbers2CSVstring(lon):
    return ",".join(str(a) for a in lon)

StringJoin = ''.join

# StringTranspose(["12345","ABCDE","abcde","54321"]) == ['1Aa5', '2Bb4', '3Cc3', '4Dd2', '5Ee1']
def StringTranspose(sList):   return [StringJoin(i) for i in zip(*sList)]

# ReplaceAll([1,"Q",77,"D"],[1,"D"],["herbivore",0.12345]) == ['herbivore', 'Q', 77, 0.12345]
def ReplaceAll(aList,a,b):
    d = dict(zip(a,b))
    return [d.get(i,i) for i in aList]


# Imports diem-specific file format: Strings of [_,0,1,2] prefixed with 'S'
def sImport(filepath):
    if isinstance(filepath,str) != True:
        raise ValueError("sImport: invalid filepath. expected str, found %s (%r)" % (type(filepath), filepath))
        #print("sImport: isinstance(filepath,str) != True: ",type(filepath),". : ",filepath)
        #return []
    with open(filepath,'r') as f:
        return [line[1:] for line in f]

# the non-prefixed equivalent of the above
def ImportString(filepath):
    with open(filepath,'r') as f:
        return f.readlines()


#def Tally(listlike):
#    return np.array(list(Counter(listlike).items()))
def HiTally(listlike):
    c = Counter(listlike)
    return np.array(list(c.most_common()))



# diem CORE helper functions
# ________________________
# csStateCount("222211100_") == [1 2 3 4]
# csStateCount(list("222211100_")) == [1 2 3 4]
def csStateCount(cs):
    """
    This function counts diem chars; input cs is char list or string; output is an np.array
    """
    ans=Counter("_012")
    ans.update(cs)
    return np.array(list(ans.values()))-1


#  pHetErrOnStateCount([1 2 3 4]) == array([0.61111111, 0.33333333, 0.1]) == [11/18,3/9,1/10]
def pHetErrOnStateCount(csCount):
    """
    This function counts diem state pState pHet pErr ; csCount is np.array return of csStateCount
    """
    sumcsCount = sum(csCount)
    if sumcsCount > 0:
        err = csCount[0]/sumcsCount
    else:
        err = np.nan
    CallTotal = sumcsCount - csCount[0]
    if CallTotal > 0:
        ans = np.array([
        sum(csCount * np.array([0,0,1,2]))/(2 * CallTotal),
        csCount[2]/CallTotal,
        err])
    else:
        ans = np.array([
        np.nan,
        np.nan,
        err])
    return ans
def pHetErrOnStringChars(cs): return pHetErrOnStateCount(csStateCount(cs))


# Here 'verbose' hook is left in intentionally: Mathematica verbose option _plots the model_
# Python version can be upgraded to do likewise
def ModelOfDiagnostic(I4, OriginalHI, Epsilon, verbose):
    if Epsilon == 0:
        return I4
    else: 
        RescaledHI = OriginalHI[~np.isnan(OriginalHI)]
        RescaledHI = (RescaledHI - np.min(RescaledHI))/(np.max(RescaledHI) - np.min(RescaledHI))
        
    SRHI = np.sort(RescaledHI)
    WarpSwitch =np.transpose( np.array([ np.mean(l2mers(SRHI),axis=1), SRHI[1:] -SRHI[:-1] ]) ) 
    WarpSwitch = sorted(WarpSwitch, key=lambda x: x[1])[-1]
    WarpSwitch = WarpSwitch[0]
    
    RescaledHIsWithNAs = MultiInsertnan(list(RescaledHI),list(np.where(OriginalHI!=OriginalHI)[0]))
    def WarpConvert(val):
        if val < WarpSwitch:
            return "0"
        elif val > WarpSwitch:
            return "2"
        else: return "_"
    WarpChars4 = list(map(WarpConvert,RescaledHIsWithNAs))
    WarpString4 = np.sum(I4[1]) * np.array(list(map(csStateCount,WarpChars4)))
    
    def pHIconvert(val):
        if val < WarpSwitch:
            return (WarpSwitch-val)/WarpSwitch
        elif val > WarpSwitch:
            return (val-WarpSwitch)/(1-WarpSwitch)
        else: return 0.5
    pHI = np.array(list(map(pHIconvert,RescaledHIsWithNAs)))
    return ((Epsilon * pHI) * WarpString4.T).T + ((1- Epsilon * pHI) * I4.T).T



# diem CORE helper functions
# CALLED IN PARALLEL
# ________________________
def GetI4ofOneCompartment2args(markerCompartmentpath, markerlabels):
    def emPolarise(markerNlabel):
        if markerNlabel[1] == "1":
            return markerNlabel[0]
        else:
            return StringReplace20(markerNlabel[0])

    Inds = StringTranspose(emPolarise(marker) for marker in zip(sImport(markerCompartmentpath), markerlabels) )
    return np.array( [csStateCount(Ind) for Ind in Inds] , dtype=np.uint32 )

def GetI4ofOneCompartment(args): return GetI4ofOneCompartment2args(*args)



# Function applying the key tight - loop function to a data compartment
def RelabelCompartmentSixArgs(markerCompartmentpath, PhiCompartment,ChosenInds,FlatLogI4,Mindex,offsetsM):
    
    # The key tight-loop function
    def PolariseAndRankMarker(origM, PhiMarker):
        Polarity = "1"
        indicesM = ReplaceAll(origM, ['_','0','1','2'], [0,1,2,3])
        revPolIndicesM = ReplaceAll(indicesM, [1,3],[3,1])
        ioM = [indicesM,offsetsM]
        offsetIndicesM    = Total(ioM)
        rioM = [revPolIndicesM,offsetsM]
        revoffsetIndicesM = Total(rioM)
        TFLI4 =    Total([FlatLogI4[index] for index in offsetIndicesM])
        revTFLI4 = Total([FlatLogI4[index] for index in revoffsetIndicesM])
#        ans = [Total(FlatLogI4[indicesM + offsetsM]),Total(FlatLogI4[revPolIndicesM + offsetsM])]
        ans = [TFLI4,revTFLI4]
        if PhiMarker == "2":
            ans = ans[::-1]
        MaxAns = max(ans)
        if ans[1] > ans[0]:
            Polarity = "2"
        PolDIsupp = [Polarity,MaxAns,MaxAns-min(ans)]
        PhiPol = [PhiMarker,Polarity]
        if   PhiPol == ["1","1"]: stringAns = ""
        elif PhiPol == ["1","2"]: stringAns = StringReplace20(origM)
        elif PhiPol == ["2","1"]: stringAns = ""
        elif PhiPol == ["2","2"]: stringAns = origM
        else: 
            print("PolariseAndRankMarker fail: PhiPol",PhiPol)
            stringAns = ""
        return [PolDIsupp,stringAns]

    AllInds = StringTranspose(sImport(markerCompartmentpath))
    MarkerStates = StringTranspose(AllInds[index] for index in ChosenInds)
    NewLabellingC = [PolariseAndRankMarker(*a) for a in zip(MarkerStates, PhiCompartment)]
    NewLabellingC = [[i] + a for i, a in enumerate(NewLabellingC)]
    ChangedStateC = [a for a in NewLabellingC if a[-1] != ""]
    DateStringAns = [time.strftime("%b %d %Y %H:%M:%S", time.localtime())]
    if len(ChangedStateC) == 0:
        DeltaIndStateCountsC = np.zeros((len(ChosenInds),4))
    else:
        IndStrings = StringTranspose([a[-1] for a in ChangedStateC])
        DeltaIndStateCountsC = [csStateCount(a) for a in IndStrings]
#        DateStringAns = [a[0] for a in ChangedStateC].insert(0,time.strftime("%b %d %Y %H:%M:%S", time.localtime()))
        DateStringAns = DateStringAns + [str(a[0]) for a in ChangedStateC]
    PolDIsupp = [a[1] for a in NewLabellingC]
#    print("DeltaIndStateCountsC[0]: ",DeltaIndStateCountsC[0])
#    print("DateStringAns[:3]: ",DateStringAns[:3])
    return (PolDIsupp, DateStringAns, DeltaIndStateCountsC)

def RelabelCompartment(args): return RelabelCompartmentSixArgs(*args)



def RePhiCompartment(PhiUpdateCompartment):
    return StringJoin(ReplaceAll(PhiUpdateCompartment, ["11","12","21","22"], ["1","2","2","1"]))



# diem results EXPORT helper functions
# CALLED IN PARALLEL
# ________________________

def ExportCompartmentPolarisation4args(compPhi,compDI,outpath,compartmentName):
    compPhiAsInt = np.array(list(compPhi)).astype(int)
    part1 = compDI[0]
    part2 = compDI[1]
    PolDIsupp = np.array(list(zip(compPhiAsInt, compDI[0], compDI[1])))
#    print("1  np.shape(PolDIsupp): ",np.shape(PolDIsupp))
#    print([type(a) for a in PolDIsupp[0]])
# auto gz does not seem to work (should be able to end file in .gz and it auto gzs)
    savepath = os.path.join(outpath, "%s.polarity.csv" % compartmentName)
    np.savetxt(savepath, PolDIsupp, fmt='%s, %s, %s')
    return savepath

def ExportCompartmentPolarisation(args): return ExportCompartmentPolarisation4args(*args)


#     PosDeltaPolarPaths = ParallelMap( ExportPosDeltaPolar, zip(STOREFULLPositionsThatChangedPolarity,list(range(len(STOREFULLPositionsThatChangedPolarity))),[outputPath] * len(STOREFULLPositionsThatChangedPolarity)) )
# EPDPfix:
def ExportPosDeltaPolar3args(compPosChangedPol,iterNumber,outpath):
# auto gz does not seem to work (should be able to end file in .gz and it auto gzs)
    savepath = os.path.join(outpath, "%d.posDeltaPolar.csv" % iterNumber)
    with open(savepath, "w") as outfile:
        outfile.write("\n".join(compPosChangedPol))
    return savepath
# END EPDPfix


def ExportPosDeltaPolar(args): return ExportPosDeltaPolar3args(*args)



# diem runtime environment setup helper functions
# ________________________

def AbsoluteTiming(f,arg):
    start_time = time.time()
    ans = f(arg)
    duration = time.time() - start_time
    return [duration,ans]

def AbsoluteTimingNoArg(f):
    start_time = time.time()
    ans = f()
    duration = time.time() - start_time
    return [duration,ans]



# diem 
# ________________________
def diem(PhiW, CompartmentNames, CompartmentPloidies, datapath, outputPath, ChosenInds, diemMaxInterations, Epsilon, verbose, processes):
 
    CompartmentChosenInds = [ChosenInds]*len(PhiW)
    I4store = []
    DIstore = []
    STOREnPositionsThatChangedPolarity = []
    STOREFULLPositionsThatChangedPolarity = []
    I4store = []
    DetectLimitCycleNchangesMatch = False
    DetectLimitCycleChangedPositionsMatch = False
    ActualLimitCycle = False,
    ExistingLCcandidate = 0
    
    SmallDataI4errorTerm = 0
    SmallDataI4errorTermStore = []
    SmallDataErrorTermGoesTo = 1
    
# Housekeeping functions
    onlyinextremeis = False
    
    vprint = print if verbose         else lambda *a: None
    xprint = print if onlyinextremeis else lambda *a: None


    def ParallelMap(f,arglist):
        with multiprocessing.Pool(processes=processes) as pool:
            ans =  pool.map(f, arglist)
            while len(pool._cache) > 0:
                sleep(0.1)
            return ans
 
    if len(CompartmentNames) == 1:
        OuterMap = map
        InnerMap = ParallelMap
    else:
        OuterMap = ParallelMap
        InnerMap = map
            
    
# Data Import functions
    def GetI4ofCompartments(markerLabelsForCompartments):
#    def GetI4ofCompartments():
        compInputs = list(zip(CompartmentDatapaths,markerLabelsForCompartments))
        return np.array(ParallelMap(GetI4ofOneCompartment,compInputs))

#        NewLabels = AbsoluteTiming(ParallelMap(RelabelCompartment,[a for a in zip(CompartmentDatapaths,Phi,CompartmentChosenInds)]))
    def RelabelCompartments():
        compInputs = list(zip(
            CompartmentDatapaths,
            Phi,
            CompartmentChosenInds,
            [FlatLogI4]*len(Phi),
            [Mindex]*len(Phi),
            [offsetsM]*len(Phi)
        ))
        compOutputs =  list(ParallelMap(RelabelCompartment,compInputs))
# and now list transpose:
        return [list(a) for a in zip(*compOutputs)]

        
    def UpdateEMstateGivenI4(locActualLimitCycle,locI4,locSmallDataI4errorTerm,locA4compartments,locA4,locFlatLogI4):
  # loc stands for UpdateEMstateGivenI4 LOCAL copy              
        vprint("UpdateEM1: Dimensions[I4]: ",locI4.shape)
        I4store.append(I4)
        if not(DetectLimitCycleChangedPositionsMatch):
            locActualLimitCycle = False
        else:
# (A==B).all()
            locActualLimitCycle = (locI4 == I4store[ExistingLCcandidate]).all()

        if not(locActualLimitCycle):
#           prepare for next iteration...
            MinInI4 = np.amin(locI4)
            locSmallDataI4errorTerm = max(0, SmallDataErrorTermGoesTo - MinInI4)
            vprint("MinInI4: ",MinInI4,". SmallDataErrorTerm:",locSmallDataI4errorTerm)
            SmallDataI4errorTermStore.append(locSmallDataI4errorTerm)
            locI4 = locI4 + locSmallDataI4errorTerm
#            locA4compartments = np.transpose(locA4compartments.T + (locSmallDataI4errorTerm * A4errorTermDistributor))
            locA4compartments = locA4compartments + (locSmallDataI4errorTerm * A4errorTermDistributor)[:,:,np.newaxis]
            print("np.shape(locA4compartments): ",np.shape(locA4compartments))
            locA4 = Total(locA4compartments)
#            print("locA4[0].shape: ",locA4[0].shape)
#            print("pHetErrOnStateCount(locA4[0]): ",pHetErrOnStateCount(locA4[0]))
            OriginalHI = np.array([pHetErrOnStateCount(a)[0] for a in locA4])
            vprint("UpdateEM2: Dimensions[I4]: ",locI4.shape)
            V4 = ModelOfDiagnostic(locI4, OriginalHI, Epsilon, verbose)
#            FlatLogI4 = N[Log[Flatten[V4/Map[Total, V4]]]]
            locFlatLogI4 = np.log(np.ndarray.flatten(np.array([a/Total(a) for a in V4])))
        else:
            vprint("GOOD HALT.")
        return (locActualLimitCycle,locI4,locSmallDataI4errorTerm,locA4compartments,locA4,locFlatLogI4)
     
    def UpdateM4withDeltaTwoArgs(M4, DeltaM4):
#        print("UpdateM4withDelta shapes:",np.shape(M4),np.shape(DeltaM4))
        TrM4 = M4.T
        TrM4[1] = TrM4[1] - DeltaM4[3] + DeltaM4[1]
        TrM4[3] = TrM4[3] - DeltaM4[1] + DeltaM4[3]
        return TrM4.T
    def UpdateM4withDelta(args): return UpdateM4withDeltaTwoArgs(*args)
    
# ACTION : Data input paths initialisation
    vprint("diem Initialisation...")
    CompartmentDatapaths = [ os.path.join(datapath, name) for name in CompartmentNames ]
# precalcs and inits
    Mindex = range(len(ChosenInds))
    offsetsM = range(0,4*len(ChosenInds)-1,4)
# ACTION : Measure initial state
    vprint("Starting big state counts..")
    I4compartments = AbsoluteTiming(GetI4ofCompartments,PhiW);
#    I4compartments = AbsoluteTimingNoArg(GetI4ofCompartments);
    vprint("I4 : ",I4compartments[0]," seconds.")
    I4compartments = I4compartments[1]
    print("np.shape(I4compartments): ",np.shape(I4compartments))
    I4errorTermDistributor = np.array([Total(a[0]) for a in I4compartments])
    print("np.shape(I4errorTermDistributor): ",np.shape(I4errorTermDistributor))
    I4errorTermDistributor = I4errorTermDistributor/Total(I4errorTermDistributor)
#   A4errorTermDistributor = CompartmentPloidies * I4errorTermDistributor
    A4errorTermDistributor = CompartmentPloidies * I4errorTermDistributor[:, np.newaxis]
    print("np.shape(A4errorTermDistributor): ",np.shape(A4errorTermDistributor))
    I4 = Total(I4compartments)
#    A4compartments = CompartmentPloidies * I4compartments
    A4compartments = CompartmentPloidies[:,:, np.newaxis] * I4compartments
    print("np.shape(A4compartments): ",np.shape(A4compartments))
    A4 = Total(A4compartments)
    print("np.shape(A4): ",np.shape(A4))
    if verbose:
        print("I4,A4 (max 5 elements shown): ")
        print(I4[:5])
        print("A4: ",A4[:5])
        
    FlatLogI4 = 0;
    ActualLimitCycle,I4,SmallDataI4errorTerm,A4compartments,A4,FlatLogI4 = UpdateEMstateGivenI4(ActualLimitCycle,I4,SmallDataI4errorTerm,A4compartments,A4,FlatLogI4)
    Phi = PhiW;
    
#   ACTION: THE EM interations:
    IterationCount = 1;
    while( not(ActualLimitCycle) and IterationCount <= diemMaxInterations):
        vprint("diem: Iteration",IterationCount)
        vprint("Starting relabelling...")
        NewLabels = AbsoluteTimingNoArg(RelabelCompartments)
        vprint("...",NewLabels[0],"seconds.")
        NewLabels = NewLabels[1]                                                
        PositionsThatChangedPolarity = [a[1:] for a in NewLabels[1]]
        xprint("iteration time stamps",IterationCount,[a[:1] for a in NewLabels[1]])
        DeltaI4compartments = NewLabels[2]
        vprint("np.shape(DeltaI4compartments): ",np.shape(DeltaI4compartments))
#        DeltaA4compartments = CompartmentPloidies * DeltaI4compartments
        DeltaA4compartments = CompartmentPloidies[:,:, np.newaxis] * DeltaI4compartments
        print("np.shape(DeltaA4compartments): ",np.shape(DeltaA4compartments))
        DeltaI4compartments = [np.array(a).T for a in DeltaI4compartments]
        DeltaA4compartments = [np.array(a).T for a in DeltaA4compartments]
        print("...post [checked valid] transpose...")
        print("np.shape(DeltaI4compartments): ",np.shape(DeltaI4compartments))
        print("np.shape(DeltaA4compartments): ",np.shape(DeltaA4compartments))
        DeltaI4 = Total(DeltaI4compartments)
        vprint("np.shape(DeltaI4): ",np.shape(DeltaI4))
        xprint("DeltaI4: ",DeltaI4)
        NewLabels = NewLabels[0]
# EPDPfix:        
        STOREFULLPositionsThatChangedPolarity.append([listOfNumbers2CSVstring(a) for a in PositionsThatChangedPolarity])
        print("***",STOREnPositionsThatChangedPolarity[-1:])
        STOREnPositionsThatChangedPolarity.append([len(a) for a in PositionsThatChangedPolarity])
        print("###",STOREnPositionsThatChangedPolarity[-1:])
        
        vprint("Lengths of NewLabels,ChangedLabels (compartment sizes, Nchanged): ")
        Remaining = np.array(list(zip([len(a) for a in NewLabels],[len(a) for a in PositionsThatChangedPolarity])))
        xprint(Remaining)
        Remaining = Total(Remaining)
        vprint("      PROGRESS (%):",100 * (1-Remaining[1]/Total(Remaining)))
        
        vprint("Extracting DeltaPhi,DI...")
        DeltaPhi = [np.array(a).T for a in NewLabels]
        DIstore.append([a[1:] for a in DeltaPhi])
    #  Spretus BODGE - this should be StringJoined BEFORE exit from constructor (WAS Map[First,[CapitalDelta][CapitalPhi]]) 
        DeltaPhi = [StringJoin(a[0]) for a in DeltaPhi]
        xprint("DeltaPhi[:5]: ",DeltaPhi[:5])
        vprint("Constructing NuPhi...")
        NuPhi = [StringTranspose(a) for a in zip(Phi,DeltaPhi)]
        xprint("NuPhi Summaries (one for each compartment): ")
#WARNINFfix        xprint(np.array([HiTally(a) for a in NuPhi]))
        vprint("starting [YES!] ParallelMap(RePhiCompartment,NuPhi)...")
        Phi = ParallelMap(RePhiCompartment,NuPhi)
#        Phi = [RePhiCompartment(a) for a in NuPhi]
        print("np.shape(Phi): ",np.shape(Phi))
        print([len(a) for a in Phi])
        print([type(a[0]) for a in Phi])
        xprint("Phi (10 chars per):")
#WARNINFfix        xprint(np.array([a[:10] for a in Phi]))
        xprint("Updated Phi Summary: ")
        xprint(np.array([HiTally(a) for a in Phi]))
        print("STOREnPos:")
        print(STOREnPositionsThatChangedPolarity)
# START Check for Actual Limit Cycle  - progressive checks to minimise computation load
#        DetectLimitCycleNchangesMatch = STOREnPositionsThatChangedPolarity[-1:] in set(STOREnPositionsThatChangedPolarity[:-1])
        ExistingLCcandidate = next((i for i, v in enumerate(STOREnPositionsThatChangedPolarity[:-1]) if v == STOREnPositionsThatChangedPolarity[-1]),-999)
        print("ExistingLCcandidate: ",ExistingLCcandidate)
        if ExistingLCcandidate >= 0:
            DetectLimitCycleNchangesMatch = True
            DetectLimitCycleChangedPositionsMatch = STOREFULLPositionsThatChangedPolarity[-1] == STOREFULLPositionsThatChangedPolarity[ExistingLCcandidate]
        else:
            DetectLimitCycleNchangesMatch = False
            DetectLimitCycleChangedPositionsMatch = False
        vprint("(DetectLimitCycleNchangesMatch,DetectLimitCycleChangedPositionsMatch,ActualLimitCycle)")
        vprint(DetectLimitCycleNchangesMatch,DetectLimitCycleChangedPositionsMatch,ActualLimitCycle)
# END Check for Actual Limit Cycle  - progressive checks to minimise computation load
#        break 
            
# UPDATE I4
        print("I4[:5]")
        print(I4[:5])
        NuI4 = UpdateM4withDelta([I4,DeltaI4])
        vprint("Comparing old and updated I4s (max 10 comparisons):")
        print("NuI4[:5]")
        print(NuI4[:5])
# UPDATE A4
#        NuA4compartments = ParallelMap(UpdateM4withDelta,zip(A4compartments,DeltaA4compartments))
        NuA4compartments = [UpdateM4withDelta(a) for a in zip(A4compartments,DeltaA4compartments)]
    
        I4 = NuI4 - SmallDataI4errorTerm
#        A4compartments = NuA4compartments - SmallDataI4errorTerm * A4errorTermDistributor;
        A4compartments = NuA4compartments - (SmallDataI4errorTerm * A4errorTermDistributor)[:,:, np.newaxis]
        A4 = Total(A4compartments)
        vprint("np.shape(A4): ",np.shape(A4))           
        
        ActualLimitCycle,I4,SmallDataI4errorTerm,A4compartments,A4,FlatLogI4 = UpdateEMstateGivenI4(ActualLimitCycle,I4,SmallDataI4errorTerm,A4compartments,A4,FlatLogI4)
        IterationCount += 1
# END While for EM interations
        """
(* Imagine we have entered a limit cycle, and we are contemplating exiting and outputing answers *)
(* Phi is updated. I4 is updated. BUT DI is NOT updated|I4' *) 
(* IFF the limit cycle is 'nothing changes' then DI' = DI and all is well...*)
(* more generally: we can take DI' from the cycle state AFTER the last occurence of the match to the 
current state in the limit cycle *)
        """
# START output to files
    vprint("Starting result exports...")
    vprint("__________________________")
    xprint("len(DIstore): ",len(DIstore))
    xprint("ExistingLCcandidate + 1: ",ExistingLCcandidate + 1)
    if ExistingLCcandidate >= 0 and ExistingLCcandidate + 1 < len(DIstore):
        vprint("Starting polarisation exports...")
        vprint("________________________________")
        CompartmentPolarisationPaths = ParallelMap( ExportCompartmentPolarisation,zip(Phi,DIstore[ExistingLCcandidate + 1],[outputPath] * len(CompartmentNames),CompartmentNames) )
        print("    Polarisation exported.")
    else:
        CompartmentPolarisationPaths = []
        print("    UNCONVERGED. NO polarisation exported. Suggest increasing diemMaxInterations.")
    if len(STOREFULLPositionsThatChangedPolarity) > 0:
        vprint("Starting PosDeltaPolar exports...")
        vprint("________________________________")
        PosDeltaPolarPaths = ParallelMap( ExportPosDeltaPolar, zip(
            STOREFULLPositionsThatChangedPolarity,
            range(len(STOREFULLPositionsThatChangedPolarity)),
            [outputPath] * len(STOREFULLPositionsThatChangedPolarity)
        ) )
        vprint("PosDeltaPolar exported.")
    else:
        PosDeltaPolarPaths = []
        vprint("    posDeltaPolar NOT exported.")   
    
    xprint("np.shape(I4): ",np.shape(I4))
    I4path = os.path.join(outputPath,"I4.tsv")
    np.savetxt(I4path, I4, delimiter = '\t')
    xprint("len(SmallDataI4errorTermStore): ",len(SmallDataI4errorTermStore))
    xprint("SmallDataI4errorTermStore: ",SmallDataI4errorTermStore)
    if len(SmallDataI4errorTermStore) > 0:
        SmallDETpath = os.path.join(outputPath,"SmallDataErrorTerm.tsv")
        np.savetxt(SmallDETpath, np.array(SmallDataI4errorTermStore), delimiter = '\t', fmt = '%10u')
    else:
        SmallDETpath = ""
    print("")
    print("diem ALL DONE.")
    alloutputpaths = [CompartmentPolarisationPaths,PosDeltaPolarPaths,[I4path,SmallDETpath]]
    return alloutputpaths

