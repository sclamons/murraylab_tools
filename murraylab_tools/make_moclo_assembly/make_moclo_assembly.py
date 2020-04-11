import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import os
import math
#import beeswarm as bs
import sys
import time
import pydna
import itertools as it
import datetime
import dnaplotlib as dpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatch
from matplotlib.patches import FancyBboxPatch
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.assembly import Assembly as pydAssembly
from Bio.Restriction import BsaI
from Bio.Restriction import BbsI
from Bio.Restriction import AarI
from Bio.Restriction import Esp3I
from copy import deepcopy as dc
import ipywidgets as widgets
from collections import defaultdict
from IPython.display import FileLink, FileLinks
import warnings
import re
def incrementString(s):
    """regular expression search! I forget exactly why this is needed"""
    m = re.search(r'\d+$', s)
    if(m):
        return s[:m.start()]+str(int(m.group())+1)
    else:
        return s+str(0)
#the following makes a few data members for handling restriction enzymes
enzymelist = [BsaI,BbsI,AarI,Esp3I]
enzymes = {str(a):a for a in enzymelist}
enlist = [str(a) for a in enzymelist]+["gibson"]
#the following defines the overhangs in our library!
ENDDICT = { \
"GGAG":"A", \
"TACT":"B", \
"AATG":"C", \
"AGGT":"D", \
"GCTT":"E", \
"CGCT":"F", \
"TGCC":"G", \
"ACTA":"H", \
"TAGA":"sc3",\
"CATTACTCGCATCCATTCTCAGGCTGTCTCGTCTCGTCTC" : "1",\
"GCTGGGAGTTCGTAGACGGAAACAAACGCAGAATCCAAGC" : "2",\
"GCACTGAAGGTCCTCAATCGCACTGGAAACATCAAGGTCG" : "3",\
"CTGACCTCCTGCCAGCAATAGTAAGACAACACGCAAAGTC" : "4",\
"GAGCCAACTCCCTTTACAACCTCACTCAAGTCCGTTAGAG" : "5",\
"CTCGTTCGCTGCCACCTAAGAATACTCTACGGTCACATAC" : "6",\
"CAAGACGCTGGCTCTGACATTTCCGCTACTGAACTACTCG" : "7",\
"CCTCGTCTCAACCAAAGCAATCAACCCATCAACCACCTGG" : "8",\
"GTTCCTTATCATCTGGCGAATCGGACCCACAAGAGCACTG" : "9",\
"CCAGGATACATAGATTACCACAACTCCGAGCCCTTCCACC" : "X",\
}
#have a dictionary of the reverse complement too
rcENDDICT = {str(Dseq(a).rc()):ENDDICT[a] for a in ENDDICT}

prevplate = None
selenzyme = "gibson" #which enzyme to assemble everything with
chewnt = 40
frags = [] #fragments in the reaction
#the following lists the components in each well, in uL. I think this is outdated
#as of 4/25/19
gga = \
[["component","volume"],
 #["buffer10x",0.4],
 #["ATP10mM",0.4],
 #["BsaI", 0.2],
 #["ligase",0.2],
 ["NEBbuffer",0.4],
 ["NEBenzyme",0.2],
 ["water",1.4],
 ["dnasln",1],
 ]
gibassy = \
[["component","volume"],
["GGAMM",1],
["dnasln",1]]
ggaPD = pd.DataFrame(gga[1:],columns=gga[0]) #this just turns it into a data frame
gibassyPD = pd.DataFrame(gibassy[1:],columns=gibassy[0])

ggaFm = 6.0
ggavecGm = 6.0
gibFm = 6.0
gibvecFm = 6.0
partsFm = ggaFm  #default is gga
vectorFm = ggavecGm
source = "384PP_AQ_BP"
ptypedict = {
            "ASSGGA04":"384PP_PLUS_AQ_BP",
            "ASSGIB01":"384LDV_PLUS_AQ_BP",
            "ASSGIB02":"384PP_AQ_BP"}
waterwell = "P1" #in your source plate, include one well that is just full of water.
#dnaPath = os.path.join(".","DNA")

#go down and look at makeEchoFile

def startText():
    print("Welcome to Moclo Assembly Helper V1")
    print("===================================")

def pickEnzyme():
    """asks the user about what kind of enzyme s/he wants to use"""
    print("Which enzyme would you like to use?")
    for el in range(len(enlist)):
                print("[{}]  {}".format(el,enlist[el]))
    print()
    userpick = int(input("type the number of your favorite! "))
    selenzyme = enlist[userpick].lower()
    print("===================================")
    return selenzyme
def findExpts(path):
    """gets a list of files/folders present in a path"""
    walkr = os.walk(path)
    dirlist = [a for a in walkr]
    expts = []
    #print(dirlist)
    #for folder in dirlist[1:]:
    folder = ['.']
    for fle in dirlist[0][2]:
        if(fle[-3:]=='csv'):
            try:
                fline = open(os.path.join(folder[0],fle),'r').readline().split(',')
                if("promoter" in fline):
                    expts+=[(os.path.join(folder[0],fle),fle[:-4])]
            except IOError:
                pass
        if(fle[-4:]=='xlsx'):
            try:
                xl_file = pd.read_excel(os.path.join(folder[0],fle),None)
                dfs = {sheet_name: xl_file.parse(sheet_name)
                          for sheet_name in xl_file.sheet_names}
                #print(dfs.keys()
                if(dfs["Sheet1"].columns[0] == "promoter"):
                    expts+=[(os.path.join(folder[0],fle),fle[:-5])]
            except (IOError,KeyError) as e:
                pass
    return sorted(expts)[::-1]

def findPartsLists(path):
    """gets a list of files/folders present in a path"""
    walkr = os.walk(path)
    dirlist = [a for a in walkr]
    #print dirlist
    expts = []
    for fle in dirlist[0][2]:
        #print fle
        if(fle[-4:]=='xlsx'):
            try:
                xl_file = pd.read_excel(os.path.join(path,fle),None)
                dfs = {sheet_name: xl_file.parse(sheet_name)
                          for sheet_name in xl_file.sheet_names}
                #print(dfs.keys()
                if("parts" in list(dfs.keys())[0]):
                    expts+=[(os.path.join(path,fle),fle[:-4])]
            except IOError:
                pass
    return sorted(expts)[::-1]

def pickPartsList():
    """user interface for picking a list of parts to use. This list must
    contain the concentration of each part as well as the 384 well location
    of each part at minimum, but better to have more stuff. Check my example
    file."""
    print("Searching for compatible parts lists...")
    pllist = findPartsLists(os.path.join(".","partslist"))
    pickedlist = ''
    if(len(pllist) <=0):
        print("could not find any parts lists :(. Make sure they are in a \
                seperate folder called 'partslist' in the same directory as this script")
    else:
        print("OK! I found")
        print()
        for el in range(len(pllist)):
            print("[{}]  {}".format(el,pllist[el][1]))
        print()
        if(len(pllist)==1):
            pickedlist = pllist[0][0]
            print("picked the only one in the list!")
        else:
            userpick = int(input("type the number of your favorite! "))
            pickedlist = pllist[userpick][0]
    openlist = pd.read_excel(pickedlist,None)
    print("===================================")
    return openlist

def pickAssembly():
    """user interface for defining assemblies to build"""
    #manual = raw_input("would you like to manually enter the parts to assemble? (y/n)")

    manual = "n"
    if(manual == "n"):
        print("searching for compatible input files...")
        time.sleep(1)
        pllist = findExpts(".")
        #print pllist
        pickedlist = ''
        if(len(pllist) <=0):
            print("could not find any assembly files")
        else:
            print("OK! I found")
            print()
            for el in range(len(pllist)):
                print("[{}]  {}".format(el,pllist[el][1]))
            print()
            if(len(pllist)==1):
                pickedlist = pllist[0][0]
                print("picked the only one in the list!")
            else:
                userpick = int(input("type the number of your favorite! "))
                pickedlist = pllist[userpick][0]
        openlist = pd.read_csv(pickedlist)
        print("===================================")
        return openlist,pickedlist
    else:
        print("sorry I haven't implemented this yet")
        pickAssembly()
    return pd.read_csv(aslist),aslist
def echoline(swell,dwell,tvol,sptype = source,spname = "Source[1]",\
                                    dpname = "Destination[1]",platebc="",partid="",partname=""):
    #if(platebc!=""):
    #    sptype = ptypedict[platebc]
    return "{},{},{},{},{},{},,,{},{},{}\n".format(spname,platebc,sptype,swell,\
                                                    partid,partname,dpname,dwell,tvol)
def echoSinglePart(partDF,partname,partfm,dwell,printstuff=True,enzyme=enzymes["BsaI"]):
    """calculates how much of a single part to put in for a number of fm."""
    try:
        pwell = partDF[partDF.part==partname].well.iloc[0]
    except IndexError:
        raise ValueError("Couldn't find the right part named '"+\
          partname+"'! Are you sure you're using the right parts list?")
        return None, None, None
    pDseq = makeDseqFromDF(partname,partDF,enzyme=enzyme)
    pconc = partDF[partDF.part==partname]["conc (nM)"]
    #concentration of said part, in the source plate
    if(len(pconc)<=0):
        #in this case we could not find the part!
        raise ValueError("Part "+part+" had an invalid concentration!"+\
                            " Are you sure you're using the right parts list?")
    pconc = pconc.iloc[0]
    pplate = partDF[partDF.part==partname]["platebc"].iloc[0]
    platet = partDF[partDF.part==partname]["platetype"].iloc[0]
    e1,e2 = echoPipet(partfm,pconc,pwell,dwell,sourceplate=pplate,sptype=platet,\
                                                    partname=partname,printstuff=printstuff)
    return e1,e2,pDseq,pplate,platet
def echoPipet(partFm,partConc,sourcewell,destwell,sourceplate=None,\
                                                partname="",sptype=None,printstuff=True):
    """does the calculation to convert femtomoles to volumes, and returns
    the finished echo line"""
    pvol = (partFm/partConc)*1000
    evol = int(pvol)
    if(evol <= 25):#im not sure what happens when the echo would round to 0.
                    #better safe than sorry and put in one droplet.
        evol = 25
    if(sourceplate==None):
        if(printstuff):
            print("===> transfer from {} to {}, {} nl".format(sourcewell,destwell,evol))
        echostring = echoline(sourcewell,destwell,evol,partname=partname)
    else:
        if(printstuff):
            print("===> transfer from {}, plate {} to {}, {} nl".format(sourcewell,sourceplate,destwell,evol))
        echostring = echoline(sourcewell,destwell,evol,spname =sourceplate,\
                            sptype= sptype,platebc = sourceplate,partname=partname)
    return echostring, evol

def makeDseqFromDF(part,partslist,col = "part",enzyme=enzymes["BsaI"]):
    """looks up the part named "part" in the column specified as col, and
    converts it into a pydna object.
    this program will check if an input sequence is a valid part.
    This involves checking a couple of things:
    1) are there only two restriction cut sites?
    2) does it have the proper overhangs?
    3) after being cut, does it produce one part with bsai sites and one part without?


    """
    pseq = partslist[partslist[col] == part].sequence.iloc[0].lower()
    pcirc = partslist[partslist[col] == part].circular.iloc[0]
    p5pover = int(partslist[partslist[col] == part]["5pend"].iloc[0])
    p3pover = int(partslist[partslist[col] == part]["3pend"].iloc[0])

    povhg = int(p5pover)
    pseqRC = str(Dseq(pseq).rc()).lower()
    if(p5pover > 0):
        pseq = pseq[p5pover:]
    elif(p5pover<0):
        pseqRC = pseqRC[:p5pover]
    if(p3pover <0):
        pseq = pseq[:p3pover]
    elif(p3pover >0):
        pseqRC = pseqRC[p5pover:]
    pDseq = Dseq(pseq,pseqRC,ovhg=povhg)
    #this defines a dsdna linear sequence
    if(pcirc):
        #this makes the sequence circular, if we have to
        pDseq = pDseq.looped()
    if(enzyme != None):

        numzymes = len(enzyme.search(pDseq,linear=not pcirc))##\
                        #len(enzyme.search(pDseq.rc(),linear=pcirc))
        if(numzymes < 2 and pcirc):
            warnings.warn("Be careful! sequence {} has only {} {} site"\
                            .format(part,numzymes,str(enzyme)))
        elif(numzymes>=2):
            try:
                testcut = pDseq.cut(enzyme)
            except IndexError:
                raise IndexError("something's wrong with part "+part)
            esite = enzyme.site.lower()
            esiterc = str(Dseq(enzyme.site).rc()).lower()
            if(numzymes > 2):
                warnings.warn("{} has {} extra {} site{}!!"\
                            .format(part,numzymes-2,str(enzyme),'s'*((numzymes-2)>1)))
            insert = []
            backbone = []
            for a in testcut:
                fpend = a.five_prime_end()
                tpend = a.three_prime_end()
                if((a.find(esite)>-1) or (a.find(esiterc)>-1)):
                    #in this case the fragment we are looking at is the 'backbone'
                    backbone+=[a]
                else:
                    #we didn't find any site sequences. this must be the insert!
                    insert+=[a]
                    if((not fpend[0]=='blunt') and \
                            (not ((fpend[1].upper() in ENDDICT) or \
                                (fpend[1].upper() in rcENDDICT)))):
                        warnings.warn("{} has non-standard overhang {}"\
                                            .format(part,fpend[1].upper()))
                    if((not tpend[0]=='blunt') and \
                            (not ((tpend[1].upper() in ENDDICT) or \
                                (tpend[1].upper() in rcENDDICT)))):
                        warnings.warn("{} has non-standard overhang {}"\
                                            .format(part,tpend[1].upper()))
            if(len(insert)==0):
                raise ValueError("{} does not produce any fragments with no cut site!".format(part))
            if(len(insert)>1):
                warnings.warn("{} produces {} fragments with no cut site".format(part,len(insert)))
            if(len(backbone)>1):
                dontwarn = False
                if(not pcirc and len(backbone)==2):
                    #in this case we started with a linear thing and so we expect it
                    #to make two 'backbones'
                    dontwarn = True
                if(not dontwarn):
                    warnings.warn("{} produces {} fragments with cut sites".format(part,len(backbone)))


    return pDseq
def bluntLeft(DSseq):
    """returns true if the left hand side of DSseq is blunt"""
    if(type(DSseq)==Dseqrecord):
        DSseq = DSseq.seq
    isblunt = (DSseq.five_prime_end()[0]=='blunt')&DSseq.linear
    return(isblunt)
def bluntRight(DSseq):
    """returns true if the right hand side of DSseq is blunt"""
    if(type(DSseq)==Dseqrecord):
        DSseq = DSseq.seq
    isblunt = (DSseq.three_prime_end()[0]=='blunt')&DSseq.linear
    return(isblunt)
def isNewDseq(newpart,partlist):
    """checks to see if newpart is contained within partlist, returns true
    if it isn't"""
    new = True
    if(type(newpart)==Dseqrecord):
        newdseqpart = newpart.seq
    #seqnewpart = str(newpart).upper()
    newcirc = newpart.circular
    #dsequid = (newpart.seq).seguid()
    #print("dsequid is "+str(dsequid))
    #dsnewpart = Dseqrecord(newpart)
    #rcnewpart = newpart.rc()
    newseguid = newdseqpart.seguid()
    #print("newseguid is "+str(newseguid))
    cseguid = None
    if(newcirc and type(newpart)==Dseqrecord):
        cseguid = newpart.cseguid()
    for part in partlist:
        if(type(part == Dseqrecord)):
            dseqpart = part.seq
        partseguid = dseqpart.seguid()

        if(newseguid==partseguid):
            new=False
            break

        #if(len(part) != len(newpart)):
            #continue
        #dspart = Dseqrecord(part)
        if(newcirc and part.circular):
            if(type(part) == Dseqrecord and cseguid != None):
                comparid = part.cseguid()
                if(comparid == cseguid):
                    new=False
                    break
            #if(seqnewpart in (str(part.seq).upper()*3)):
            #    new=False
            #    break
            #elif(seqnewpart in (str(part.seq.rc()).upper()*3)):
            #    new=False
            #    break
        #elif(part == newpart or part == rcnewpart):
            #new=False
            #break
    return new
def allCombDseq(partslist,resultlist = []):
    '''recursively finds all possible paths through the partslist'''
    if(len(partslist)==1):
        #if there's only one part, then "all possible paths" is only one
        return partslist
    else:
        #result is the final output
        result = []
        for p in range(len(partslist)):
            newplist = dc(partslist)
            #basically the idea is to take the first part,
            #and stick it to the front of every other possible assembly
            part = newplist.pop(p)
            #this is the recursive part
            prevresult = allCombDseq(newplist)
            partstoadd = []
            freezult = dc(result)
            #for z in prevresult:

            for b in prevresult:
                #maybe some of the other assemblies
                #we came up with in the recursive step
                #are the same as assemblies we will come up
                #with in this step. For that reason we may
                #want to cull them by not adding them
                #to the "parts to add" list
                if(isNewDseq(b,freezult)):
                    partstoadd+=[b]
                #try to join the given part to everything else
                if((not bluntRight(part)) and (not bluntLeft(b)) and part.linear and b.linear):
                    #this means we don't allow blunt ligations! We also don't allow
                    #ligations between a linear and a circular part. Makes sense right?
                    #since that would never work anyway
                    newpart = None
                    try:
                        #maybe we should try flipping one of these?
                        newpart= part+b

                    except TypeError:
                        #this happens if the parts don't have the right sticky ends.
                        #we can also try rotating 'part' around
                        pass
                    try:
                        #part b is not blunt on the left so this is OK,
                        #since blunt and not-blunt won't ligate
                        newpart = part.rc()+b
                    except TypeError:
                        pass
                    if(newpart == None):
                        #if the part is still None then it won't ligate forwards
                        #or backwards. Skip!
                        continue
                    try:
                        if((not bluntRight(newpart)) and (not bluntLeft(newpart))):
                            #given that the part assembled, can it be circularized?
                            newpart = newpart.looped()
                            #this thing will return TypeError if it can't be
                            #looped
                    except TypeError:
                        #this happens if the part can't be circularized
                        pass
                    if(isNewDseq(newpart,result)):
                        #this checks if the sequence we just made
                        #already exists. this can happen for example if we
                        #make the same circular assembly but starting from
                        #a different spot around the circle
                        result+=[newpart]
            result+=partstoadd
        return result

def pushDict(Dic,key,value):
    """adds a value to a dictionary, whether it has a key or not"""
    try:
        pval = Dic[key]
    except KeyError:
        if(type(value)==list or type(value)==tuple):
            value = tuple(value)
            pval = ()
        elif(type(value)==str):
            pval = ""
        elif(type(value)==int):
            pval = 0
        elif(type(value)==float):
            pval = 0.0
    Dic[key] =pval + value
def findFilesDict(path=".",teststr = "promoter"):
    """gets a list of files/folders present in a path"""
    walkr = os.walk(path)
    dirlist = [a for a in walkr]
    expts = {}
    #print(dirlist)
    #for folder in dirlist[1:]:
    folder = [path]
    #print(dirlist)
    for fle in dirlist[0][2]:
        if(fle[-3:]=='csv'):
            try:
                #print('{}\\{}'.format(folder[0],fle))
                fline = open(os.path.join(folder[0],fle),'r').readline().split(',')
                if(teststr in fline):
                    expts[fle[:-4]]=os.path.join(folder[0],fle)
            except IOError:
                pass
        if(fle[-4:]=='xlsx'):
            try:
                xl_file = pd.read_excel(os.path.join(folder[0],fle))
                #dfs = {sheet_name: xl_file.parse(sheet_name)
                #          for sheet_name in xl_file.sheet_names}
                #print(dfs.keys()
                #print(xl_file.columns)
                if(teststr in xl_file.columns):
                    #print("found")
                    expts[fle[:-5]]=os.path.join(folder[0],fle)
            except (IOError,KeyError) as e:
                pass
    return expts
def findPartsListsDict(path,teststr = "parts_1"):
    """gets a list of files/folders present in a path"""
    walkr = os.walk(path)
    dirlist = [a for a in walkr]
    #print(dirlist[0][2])
    expts = {}
    for fle in dirlist[0][2]:
        #print fle
        if((fle[-4:]=='xlsx') or (fle[-4:]=='xlsm')):
            try:
                dfs = pd.read_excel(os.path.join(path,fle),None)
                #dfs = {sheet_name: xl_file.parse(sheet_name)
                #          for sheet_name in xl_file.sheet_names}
                #print(dfs)
                #print(dfs.keys())
                if(teststr in list(dfs.keys())[0]):
                    expts[fle[:-5]] = os.path.join(path,fle)
            except IOError:
                pass
    return expts

def findDNAPaths(startNode,nodeDict,edgeDict):
    """given a start, a dictionary of nodes, and a dictionary of edges,
    find all complete paths for a DNA molecule
    Complete is defined as: producing a molecule with all blunt edges,
    or producing a circular molecule."""
    #we assemble the DNA sequences from left to right.
    nnode = dc(nodeDict)
    noderight = nnode[startNode][1] #the right-hand overhang of the node in question.
    del nnode[startNode]
    destinations = edgeDict[noderight] #this could contain only one entry, the starting node
    seqs = [] #haven't found any complete paths yet
    nopaths = True
    candidateSeqs = []
    if(noderight != "blunt"): #blunt cannot go on
        for destination in destinations:
            #go through the list of destinations and see if we can go forward
            if(destination[1]==0): #this node links to something else
                if(destination[0] in nnode): #we havent visited it yet
                    nopaths = False

                    newpaths = findDNAPaths(destination[0],nnode,edgeDict) #find all paths from there!
                    for path in newpaths:
                        candidateSeqs+=[[startNode]+path]
    if(nopaths): #if we dont find any paths, call it good
        candidateSeqs+=[[startNode]]
    #print("canseqs is {}".format(candidateSeqs))
    return candidateSeqs
def getOverhang(Dnaseq,side="left"):
    """extracts the overhang in the DNA sequence, either on the left or right sides.
    If the dna sequence is blunt, then the returned overhang is called 'blunt'"""

def appendPart(part,pind,edgeDict,nodeDict):
    """this function appends a part to a dictionary of
    edges (overhangs), and nodes(middle sequence) for running DPallcombDseq.
    part is a DseqRecord of a DNA part that's been cut by an enzyme.
    pind is the index of that part in the parts list
    edgedict is a dictionary of edges that says which nodes they are connected
    to.
    nodedict is a dictionary of nodes that says which edges they have."""
    Lend = ""
    Rend = ""
    Ltype,Lseq = part.five_prime_end()
    Rtype,Rseq = part.three_prime_end()
    if(Ltype == "blunt"):
        Lend = "blunt"
        #if the end is blunt append nothing
        edgeDict[Lend].append([pind,0])
        #pushDict(edgeDict,Lend,((pind,0),))
    else:
        if(Ltype == "3'"):
            #if we have a 3' overhang, then add that sequence
            Lend = str(Dseq(Lseq).rc()).lower()
        else:
            #otherwise, it must be a 5' overhang since we handled the
            #blunt condition above.
            Lend = str(Lseq).lower()
        edgeDict[Lend].append([pind,0])
    if(Rtype == "blunt"):
        #same thing for the right side
        Rend = "blunt"
        edgeDict[Rend].append([pind,1])
    else:
        if(Rtype == "5'"):
            Rend = str(Dseq(Rseq).rc()).lower()
        else:
            Rend = str(Rseq).lower()
        edgeDict[Rend].append([pind,1])
    nodeDict[pind] = (Lend,Rend)
def annotateScar(part, end='3prime'):
    plen = len(part)
    if(end=='3prime'):
        ovhg = part.seq.three_prime_end()
        loc1 = plen-len(ovhg[1])
        loc2 = plen
    else:
        ovhg = part.seq.five_prime_end()
        loc1 = 0
        loc2 = len(ovhg[1])
    oseq = str(ovhg[1]).upper()
    scarname = "?"
    floc = int(loc1)
    sloc = int(loc2)
    dir = 1
    #scardir = "fwd"
    if((oseq in ENDDICT.keys()) or (oseq in rcENDDICT.keys())):
        #either direction for now...
        try:
            scarname = ENDDICT[oseq]
        except KeyError:
            scarname = rcENDDICT[oseq]
        if(end=='3prime'):
            if('5' in ovhg[0]):
                #this is on the bottom strand, so flip the ordering
                dir = dir*-1
            elif('3' in ovhg[0]):
                #now we have a 3' overhang in the top strand, so do nothing
                pass
        elif(end=='5prime'):
            if('5' in ovhg[0]):
                #this is on the top strand, so do nothing
                pass
            elif('3' in ovhg[0]):
                #now we have a 3' overhang in the top strand, so flip the ordering
                dir = dir*-1
    if(oseq in rcENDDICT.keys()):
        #so if we found the reverse complement in fact, then reverse everything
        #again
        dir = dir*-1
    if(dir==-1):
        floc = int(loc2)
        sloc = int(loc1)
    #oseq = str(Dseq(oseq).rc())
    part.add_feature(floc,sloc,label=scarname,type="Scar")
def DPallCombDseq(partslist):
    '''Finds all paths through the partsist using a graph type of approach.
    First a graph is constructed from all possible overhang interactions,
    then the program makes paths from every part to a logical conclusion
    in the graph, then it backtracks and actually assembles the DNA.'''
    #actually, we need to produce a graph which describes the parts FIRST
    #then, starting from any part, traverse the graph in every possible path and store
    #the paths which are "valid" i.e., produce blunt ended or circular products.
    edgeDict = defaultdict(lambda : []) #dictionary of all edges in the partslist!
    nodeDict = {}#defaultdict(lambda : [])
    partDict = {}#defaultdict(lambda : [])
    pind = 0
    import time
    rcpartslist = []
    number_of_parts = len(partslist)
    for part in partslist:
        #this next part appends the part to the list of nodes and edges
        appendPart(part,pind,edgeDict,nodeDict)
        appendPart(part.rc(),pind+number_of_parts,edgeDict,nodeDict)
        rcpartslist+=[part.rc()]
        pind+=1
    partslist+=rcpartslist
    paths = []
    for pind in list(nodeDict.keys()):
        #find good paths through the graph starting from every part
        paths += findDNAPaths(pind,nodeDict,edgeDict)
    goodpaths = []
    part1time = 0
    part2time = 0
    for path in paths:
        #here we are looking at the first and last parts
        #to see if they are blunt
        fpart = path[0]
        rpart = path[-1]
        npart = False
        accpart = Dseqrecord(partslist[fpart])
        if(nodeDict[fpart][0]=="blunt" and nodeDict[rpart][1]=="blunt"):
            #this means we have a blunt ended path! good
            npart = True
            plen = len(accpart)
            #accpart.add_feature(0,3,label="?",type="scar")
            #accpart.add_feature(plen-4,plen,label="?",type="scar")
            for pind in path[1:]:
                #this traces back the path
                #we want to add features as we go representing the cloning
                #scars. These scars could be gibson or golden gate in nature
                #SCARANNOT
                '''
                ovhg = accpart.seq.three_prime_end()
                oseq = ovhg[1]
                plen = len(accpart)
                if("5" in ovhg[0]):
                    #ideally we take note of what type of overhang it is
                    #but for now i'll just take the top strand sequence
                    oseq = str(Dseq(oseq).rc())
                accpart.add_feature(plen-len(oseq),plen,label="?",type="scar")
                #/scarannot'''
                annotateScar(accpart)
                accpart+=partslist[pind]


        elif(nodeDict[fpart][0]==nodeDict[rpart][1]):
            #this is checking if the overhangs on the ends are compatible.
            #if true, then create a circular piece of DNA!
            npart = True
            #this means we have a circular part! also good!
            #accpart = partslist[fpart]
            for pind in path[1:]:
                #SCARANNOT
                '''
                ovhg = accpart.seq.three_prime_end()
                oseq = ovhg[1]
                plen = len(accpart)
                if("5" in ovhg[0]):
                    #ideally we take note of what type of overhang it is
                    #but for now i'll just take the top strand sequence
                    oseq = str(Dseq(oseq).rc())
                accpart.add_feature(plen-len(oseq),plen,label="?",type="scar")
                #/scarannot'''
                annotateScar(accpart)
                accpart+=partslist[pind]
            #SCARANNOT
            '''
            ovhg = accpart.seq.three_prime_end()
            oseq = ovhg[1]
            plen = len(accpart)
            if("5" in ovhg[0]):
                #ideally we take note of what type of overhang it is
                #but for now i'll just take the top strand sequence
                oseq = str(Dseq(oseq).rc())
            accpart.add_feature(plen-len(oseq),plen,label="?",type="scar")
            #/scarannot'''
            annotateScar(accpart)
            accpart=accpart.looped()
        if(npart):
            #this checks if the part we think is good already exists
            #in the list
            if(isNewDseq(accpart,goodpaths)):
                goodpaths+=[accpart]
        #part2time+=time.time()-stime
    #dtime = time.time()-stime
    #stime = time.time()
    #print("done tracing back paths, took "+str(dtime))
    #print("first half took " + str(part1time))
    #print("second half took " + str(part2time))
    return goodpaths
def chewback(seqtochew,chewamt,end="fiveprime"):
    """chews back the amount mentioned, from the end mentioned."""
    wat = seqtochew.watson
    cri = seqtochew.crick

    if(len(seqtochew) > chewamt*2+1):
        if(end=="fiveprime"):
            cwat = wat[chewamt:]
            ccri = cri[chewamt:]

        else:
            cwat = wat[:-chewamt]
            ccri = cri[:-chewamt]
        newseq = Dseq(cwat,ccri,ovhg = chewamt)
        return newseq
    else:
        return None

def makeEchoFile(parts,aslist,gga=ggaPD,partsFm=partsFm,source=source,\
            output = "output.csv",selenzyme=selenzyme,fname="recentassembly",\
            protocolsDF=None,sepfiles=True,sepfilename="outputLDV.csv",\
            printstuff=True,progbar=None,mypath=".",annotateDF=None):
    """makes an echo csv using the given list of assemblies and source plate of
    parts..
    inputs:
        parts: dataframe of what's in the source plate
        aslist: dataframe of what we need to assemble
        gga: a short dictionary indicating what volume of all the components
            go into the reaction mix
        partsFm: how many femtomoles of each part to use
        source: the name of the source plate. like "384PP_AQ_BP or something
        output: the name of the output file
        selenzyme: the enzyme we are going to use for assembly. everything
            is assembled with the same enzyme! actually this does nothing because
            the enzyme is taken from the aslist thing anyway
        fname: this is the name of the folder to save the successfully assembled
            dna files into
        protocolsDF: a dataframe containing a descriptor for different possible
            protocols. For instance it would say how much DNA volume and
            concentration we need for GGA or gibson."""

    #this is the boilerplate columns list
    dnaPath = os.path.join(mypath,"DNA")
    outfile = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,\
    Sample ID,Sample Name,Sample Group,Sample Comment,Destination Plate Name,\
    Destination Well,Transfer Volume\n"
    f1init = len(outfile)
    outfile2 = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,\
    Sample ID,Sample Name,Sample Group,Sample Comment,Destination Plate Name,\
    Destination Well,Transfer Volume\n"
    f2init = len(outfile2)
    #this iterates through rows in the assembly list file. Each row
    #defines an assembly, with the columns representing what parts go in.
    #this may not be ideal but it's fairly human readable and we only do
    #four parts + vector for each assembly.
    _,fname = os.path.split(fname)
    if("." in fname):
        fname = fname[:fname.index(".")]

    #the following is for making a spreadsheet style sequence list for
    #performing further assemblies
    prodSeqSpread = "well,part,description,type,left,right,conc (nM),date,numvalue,sequence,circular,5pend,3pend,length\n"
    prevplate = None
    prevtype = None
    maxprog = float(len(aslist))

    for assnum in range(len(aslist)):
        #this goes row by row
        if(progbar != None):
            progbar.value=float(assnum+1)/maxprog
        assembly = aslist[assnum:assnum+1] #cuts out one row of dataframe
        dwell = assembly.targwell[assembly.targwell.index[0]] #well where assembly will happen

        #print("pick enzyme")
        #print(assembly)
        enzyme=None
        #if we are doing Gibson assembly, then the restriction enzyme is undefined
        try:
            selenzyme = assembly.enzyme[assembly.enzyme.index[0]]
            #if the user forgot to define an enzyme assume it is BsaI. That's the most common one we use
        except KeyError:
            selenzyme = "BsaI"
        if(protocolsDF!=None):
            cprt_temp = "gga"
            if(selenzyme == "gibson"):
                cprt_temp = "gibson"
            #iloc[0] is used in case there are multiple parts with the same
            #name. Only the first one is used in that case.
            curprot = {"dnasln": protocolsDF[(protocolsDF.protocol==cprt_temp)&\
                            (protocolsDF.component == "dnasln")].amount.iloc[0]}
            partsFm = curprot[curprot.component==partfm].amount.iloc[0]
            vectorFm = curprot[curprot.component==vectorfm].amount.iloc[0]
        else:
            curprot = ggaPD
            partsFm = ggaFm
            vectorFm = ggavecGm
            if(selenzyme == "gibson"):
                #for gibson assembly the protocol is different
                curprot = gibassyPD
                partsFm = gibFm
                vectorFm = gibvecFm
        water = float(curprot[curprot.component=="dnasln"].volume)*1000 #total amount of water, to start with
        if(printstuff):
            print("assembling with "+selenzyme)
        aind = assembly.index[0] #necessary for dataframes probably because I'm dumb
        frags = []
        if(not selenzyme == "gibson"):
            enzyme = enzymes[selenzyme]
            esite = enzyme.site.lower()
            esiterc = str(Dseq(enzyme.site).rc()).lower()
        for col in assembly:
            if(col=="targwell"):#since every row is terminated by the "target well",
                                #we'll take this opportunity to put in the water
                if(int(water) <25):
                    #echo gets mad if you tell it to pipet significantly less than 25 nl
                    water = 25
                ewat = int(water) #the echo automatically rounds to the nearest 25,
                                #so it's not really necessary to round here.
                #dsrfrags = [Dseqrecord(a) for a in frags]
                #x = pydAssembly(dsrfrags,limit = 4)
                #print(frags)
                #print(len(frags))
                allprod= []
                nefrags = []
                cutfrags = []
                if(selenzyme != "gibson"):
                    enzyme = enzymes[selenzyme]
                for frag in frags:
                    if(selenzyme == "gibson"):
                        if(len(frag)>chewnt*2+1):
                            nefrags += [chewback(frag,chewnt)]
                        else:
                            raise ValueError("part with sequence "+frag+" is too "+\
                                            "short for gibson! (<= 80 nt)")
                    else:
                        newpcs = frag.cut(enzyme)
                        if(len(newpcs) == 0):
                            newpcs+=[frag]
                        for pcs in newpcs:
                            if(pcs.find(esite)+pcs.find(esiterc)==-2):
                                nefrags+=[pcs]
                allprod = DPallCombDseq(nefrags)
                if(printstuff):
                    print("found {} possible products".format(len(allprod)))
                goodprod = []
                newpath = os.path.join(dnaPath,fname)
                if(printstuff):
                    print("saving in folder {}".format(newpath))
                Cname = ""
                try:
                    #this part gathers the "name" column to create the output sequence
                    Cname = assembly.name[assembly.name.index[0]]
                except KeyError:
                    Cname = ""
                if(Cname == "" or str(Cname) == "nan"):
                    Cname = "well"+dwell
                if(printstuff):
                    print("Parts in construct {}".format(Cname))
                if not os.path.exists(newpath):
                    if(printstuff):
                        print("made dirs!")
                    os.makedirs(newpath)

                num = 0
                for prod in allprod:
                    Cnamenum = Cname
                    #filename = Cname+".gbk"
                    if(len(allprod) > 1):
                        #filename = Cname+"_"+str(num)+".gbk"
                        #wout = open(os.path.join(newpath,filename),"w")
                        Cnamenum = Cname+"_"+str(num)
                    else:
                        pass
                        #wout = open(os.path.join(newpath,filename),"w")
                    if((bluntLeft(prod) and bluntRight(prod)) or (prod.circular)):
                        num+=1
                        goodprod+=[prod]
                        #topo = ["linear","circular"][int(prod.circular)]
                        booltopo = ["FALSE","TRUE"][int(prod.circular)]
                        #wout.write("\r\n>Construct"+str(num)+"_"+topo)
                        un_prod = "_".join(Cnamenum.split())
                        #wout.write("LOCUS       {}                {} bp ds-DNA     {} SYN 01-JAN-0001\n".format(un_prod,len(prod),topo))
                        #wout.write("ORIGIN\n")
                        #wout.write(str(prod)+"\n//")
                        now = datetime.datetime.now()
                        nowdate = "{}/{}/{}".format(now.month,now.day,now.year)
                        prod.name = Cnamenum
                        plt.figure(figsize=(8,1))
                        ax = plt.gca()
                        drawConstruct(ax,prod,annotateDF=annotateDF)
                        plt.show()
                        prod.write(os.path.join(newpath,Cnamenum+".gbk"))
                        prodSeqSpread += "{},{},assembled with {},,,,30,{},,{},{},{},{},{}\n".format(\
                                        dwell,un_prod,          selenzyme,nowdate,prod.seq,booltopo,0,0,len(prod))
                    #wout.close()
                assembend = ["y","ies"][int(len(goodprod)>1)]
                if(printstuff):
                    print("Detected {} possible assembl{}".format(len(goodprod),assembend))
                frags = []
                if(water <=0):
                    print("WARNING!!!! water <=0 in well {}".format(dwell))
                else:
                    #print("water from {} to {}, {} nl".format(waterwell,dwell,ewat))
                    if(prevplate == None):
                        #print("normalwater")
                        #im not convinced this ever gets triggered
                        #but just in case, i guess we can find the first water well
                        waterrows=parts[parts.part=="water"]
                        if(len(waterrows)==0):
                            raise KeyError("no water wells indicated!")
                        #print(waterrows)
                        waterrow = waterrows.iloc[0]
                        waterwell = waterrow.well
                        platetype= waterrow.platetype
                        curplatebc = waterrow.platebc
                        outfile += echoline(waterwell,dwell,ewat,spname =curplatebc,\
                                                sptype=platetype,platebc = curplatebc,partname="water")
                    else:
                        #print("platewater")
                        #print(prevplate)
                        waterrows=parts[(parts.part=="water") & (parts.platebc==prevplate)]
                        if(len(waterrows)==0):
                            raise KeyError("no water wells indicated!")
                        #print(waterrows)
                        waterrow = waterrows.iloc[0]
                        waterwell = waterrow.well
                        watline = echoline(waterwell,dwell,ewat,spname =prevplate,\
                                                sptype=prevtype,platebc = prevplate,partname="water")
                        if("LDV" in prevtype):
                            outfile2+=watline
                        else:
                            outfile += watline
                    #add water to the well!
                if(printstuff):
                    print("")
            elif(col in ["comment","enzyme","name"]):#skip this column!
                pass
            else:
                #this is the part name from the "assembly" file
                part = assembly[col][aind]

                if(str(part) == 'nan'):
                    #this means we skip this part, because the name is empty
                    if(printstuff):
                        print("skip one!")
                else:
                    #shouldnt need to define "part" again??
                    #part = assembly[col][aind]
                    #this is the name of the part!
                    #parts[parts.part==assembly[col][aind]].well.iloc[0]
                    evol = 0
                    if(':' in str(part)):
                        #this means we have multiple parts to mix!
                        subparts = part.split(':')
                        t_partsFm = partsFm/len(subparts)
                        t_vecFm = vectorFm/len(subparts)
                        for subpart in subparts:
                            useFm = t_partsFm
                            if(col == "vector"):
                                #use the vector at lower concentration!!
                                useFm = t_vecFm
                            e1,e2,pDseq,prevplate,prevtype = echoSinglePart(parts,\
                                    subpart,useFm,dwell,printstuff=printstuff,enzyme=enzyme)
                            frags+=[pDseq]
                            evol += e2
                            if(sepfiles):
                                if("LDV" in e1):
                                    outfile2+=e1
                                else:
                                    outfile+= e1
                            else:
                                outfile+= e1


                    else:
                        useFm = partsFm
                        if(col == "vector"):
                            #use the vector at lower concentration!!
                            useFm = vectorFm
                        e1,e2,pDseq,prevplate,prevtype = echoSinglePart(parts,\
                                    part,useFm,dwell,printstuff=printstuff,enzyme=enzyme)
                        frags+=[pDseq]
                        evol += e2
                        if(sepfiles):
                            if("LDV" in e1):
                                outfile2+=e1
                            else:
                                outfile+= e1
                        else:
                            outfile+= e1
                    water=water-evol
    pspread = open(os.path.join(newpath,fname+".csv"),"w")
    pspread.write(prodSeqSpread)
    pspread.close()
    seqdispDF = pd.read_csv(os.path.join(newpath,fname+".csv"),usecols=["well","part","circular","length"])
    display(seqdispDF)
    display(FileLink(os.path.join(newpath,fname+".csv")))
    if(len(outfile)>f1init):
        ofle = open(output,"w")
        ofle.write(outfile)
        ofle.close()
        display(FileLink(output))
    if(sepfiles and (len(outfile2) > f2init)):
        if(printstuff):
            print("wrote LDV steps in {}".format(sepfilename))
        ofle2 = open(sepfilename,"w")
        ofle2.write(outfile2)
        ofle2.close()
        display(FileLink(sepfilename))
outitems = []


class assemblyFileMaker():
    def __init__(self,mypath=".",partsdf = None):
        self.p = partsdf
        self.holdup=False
        self.ddlay = widgets.Layout(width='75px',height='30px')
        self.eblay = widgets.Layout(width='50px',height='30px')
        self.lsblay = widgets.Layout(width='140px',height='30px')
        self.sblay = widgets.Layout(width='100px',height='30px')
        self.rsblay = widgets.Layout(width='60px',height='30px')
        self.Vboxlay = widgets.Layout(width='130px',height='67px')
        self.textlay = widgets.Layout(width='200px',height='30px')
        self.PlateLetters="ABCDEFGHIJKLMNOP"
        self.PlateNumbers=(1,2,3,4,5,6,7,8,9,10,11,12,\
                                13,14,15,16,17,18,19,20,21,22,23,24)
        self.PlateRowsCols=(16,24)
        self.mypath = mypath
        if(type(self.p)==pd.DataFrame):
            self.parts={"google doc":"google doc"}
        else:
            self.parts = findPartsListsDict(os.path.join(self.mypath,"partslist"))

        #txtdisabl = False
        assemblies = []
        oplist = findFilesDict(os.path.join(mypath,"assemblies"))
        #parts = findPartsListsDict(os.path.join(mypath,"partslist"))

        self.loadFIleList = widgets.Dropdown(
            options=oplist,
            #value=2,
            layout=self.lsblay,
            description='',
        )
        self.loadbut = widgets.Button(
            description='Load',
            disabled=False,
            button_style='', # 'success', 'info', 'warning', 'danger' or ''
            layout=self.rsblay,
            tooltip='Click to load an existing file',
        )




        self.listEverything = widgets.Checkbox(
            value=False,
            description='List all parts',
            disabled=False
        )
        self.fname1 = widgets.Text(
            value="untitled",
            placeholder = "type something",
            description='Assembly File Name:',
            layout=self.textlay,
            disabled=False
        )
        self.DestWell = widgets.Text(
            value="A1",
            placeholder = "type something",
            description='Dest Well:',
            layout=self.Vboxlay,
            disabled=True
        )
        self.AddCols = widgets.IntText(
            value=0,
            placeholder = "type something",
            description='Extra Cols:',
            layout=self.Vboxlay,
            #disabled=True
        )
        self.drop2 = widgets.Dropdown(
            options=self.parts,
            width=100,
            #value=2,
            description='parts list:',
            layout=self.textlay,
        )
        #print(self.drop2.style.keys)
        self.but = widgets.Button(
            description='New...',
            disabled=False,
            button_style='', # 'success', 'info', 'warning', 'danger' or ''
            layout=self.sblay,
            tooltip='Click to start adding assemblies',

            #icon='check'
        )
        self.finbut = widgets.Button(
            description='Save!',
            disabled=True,
            button_style='warning',#, 'danger' or ''
            layout=self.sblay,
            tooltip='Finish and Save',

            #icon='check'
        )

        self.but.on_click(self.on_button_clicked)
        self.finbut.on_click(self.finishAndSave)
        self.loadbut.on_click(self.loadFile_clicked)
        self.listEverything.observe(self.on_listEverything_changed,names='value')
        self.cbox = widgets.HBox([
                    widgets.VBox([self.fname1,widgets.HBox([self.loadFIleList,self.loadbut]),self.listEverything]),\
                    widgets.VBox([self.drop2,widgets.HBox([self.DestWell,self.AddCols])]),\
                    widgets.VBox([self.but,self.finbut],layout=self.Vboxlay)])
        display(self.cbox)
    def add_row(self,b):
        thisrow = int(b.tooltip[4:])
        self.addWidgetRow(labonly=False,copyrow=thisrow)
        outcols = [widgets.VBox(a) for a in self.outitems ]
        self.bigSheet.children=outcols
        #b.disabled=True
        #print(b)
    def remove_row(self,b):
        thisrow = int(b.tooltip[4:])
        #outcolnum=0
        cleared = False
        for colnum in list(range(len(self.outitems))[:-3])\
                                                    +[len(self.outitems)-2]:
            pvalue = self.outitems[colnum][thisrow].value
            if(pvalue != ""):
                cleared = True
            self.outitems[colnum][thisrow].value = ""
        if(cleared):
            return

        for colnum in range(len(self.outitems)):
            self.outitems[colnum]=self.outitems[colnum][:thisrow]+\
                        self.outitems[colnum][thisrow+1:]

            #outcolnum +=1
        newbutcol = []
        newrow = 0
        for a in self.outitems[-1]:
            #print(a)
            try:
                a.children[0].tooltip = "row "+str(newrow)
                a.children[1].tooltip = "row "+str(newrow)
                if(len(self.outitems[0])<=2):
                    a.children[1].disabled=True
                else:
                    a.children[1].disabled=False
            except AttributeError:
                pass
            newrow +=1
        outcols = [widgets.VBox(a) for a in self.outitems ]
        self.bigSheet.children=outcols
        #print(b)
    def generateOptionsList(self,df,colname,prevval=None,listmode=0):
        """come up with a list of options given a column name. This contains
        a ton of specific code"""
        oplist = []
        if(listmode == 1 and colname != "enzyme"):
            oplist = sorted(list(df.part))+[""]
        else:
            if("vector" in colname):
                oplist = sorted(list(df[(df.type=="UNS")|\
                                        (df.type=="vector")].part))+[""]
            elif(colname=="enzyme"):
                oplist =enlist
                if(prevval == ""):
                    prevval = enlist[0]
            else:
                oplist = sorted(list(df[df.type==colname].part))+[""]
        if(not (prevval in oplist)):
            oplist+=[prevval]
        return oplist,prevval
    def on_listEverything_changed(self,change):
        """this triggers when you change the value of "listEverything".
        Here we want to change the values in the drop down to correspond to
        either
        (a) surrounding parts or
        (b) the appropriate category
        """
        self.updatePartOptions(None)
        """
        typewewant = type(widgets.Dropdown())
        #this means we checked the box. Now change drop box's options
        for col in self.outitems:
            for item in col:
                if(type(item)==typewewant):
                    oplist,pval = self.generateOptionsList(self.p,\
                                        col[0].value,item.value,change['new'])
                    item.options=oplist
                    item.value=pval
        #"""
    def loadFile_clicked(self,b):
        """loads a file from memory, instead of making a brand new one!"""
        self.on_button_clicked(b,loadFile=self.loadFIleList.value)
    def on_button_clicked(self,b,loadFile=None):
        """start making the assembly! THis part loads the first row of parts
        drop downs and populates them with options!"""
        #txtdisabl = True
        b.disabled=True
        self.but.disabled = True
        self.drop2.disabled = True
        self.finbut.disabled = False
        self.DestWell.disabled = False
        self.AddCols.disabled = True
        self.loadFIleList.disabled=True
        self.loadbut.disabled=True
        if(loadFile!=None):
            #this should read the file
            self.fname1.value=os.path.splitext(os.path.split(loadFile)[1])[0]
            ftoload = pd.read_csv(loadFile).fillna('')
            try:
                ftoload = ftoload.drop('comment',axis=1)
            except (ValueError,KeyError) as e:
                #if this happens then 'comment' was already not there. great!
                pass

            self.AddCols.value=len(ftoload.columns)-9

        if(not(type(self.p)==pd.DataFrame)):
            dfs = pd.read_excel(self.drop2.value,None)
            sheetlist = list(dfs.keys())
            self.p = pd.DataFrame.append(dfs["parts_1"],dfs["Gibson"])
        self.collabels = ["vector1","promoter","UTR","CDS","Terminator","vector2","enzyme","name",""]
        if(self.AddCols.value>0):
            newclabeld = self.collabels
            for x in range(self.AddCols.value):
                newclabeld=newclabeld[:-4]+["newcol"+str(x+1)]+newclabeld[-4:]
            self.collabels = newclabeld
        self.outitems = []
        self.addWidgetRow(labonly=True)
        if(loadFile==None):
            self.addWidgetRow(labonly=False)
        else:
            #print(loadFile)

            findex = ftoload.index
            first = True
            for findex in ftoload.index:
                dfrow = ftoload.iloc[findex]
                currow = list(dfrow)
                if(first):
                    self.DestWell.value=dfrow.targwell
                    #extracols =
                    #startpos =
                    first=False
                currow = list(dfrow.drop(['targwell','name','enzyme']))\
                                +[dfrow.enzyme]+[dfrow["name"]]
                self.addWidgetRow(labonly=False,copyrow=currow)
            #self.updatePartOptions()
            #readindex = ftoload.index()
        outcols = [widgets.VBox(a) for a in self.outitems ]
        self.bigSheet=widgets.HBox(outcols)
        display(self.bigSheet)
    def updatePartOptions(self,b=None):
        """update the options available to each drop down, according to what
        values are chosen in the other drop downs. For example, only allow
        parts which are compatible"""
        if(self.holdup):
            return
        self.holdup=True
        maxcols = len(self.outitems)-3
        for colnum in range(maxcols):
            for itemnum in range(len(self.outitems[colnum]))[1:]:
                curitem = self.outitems[colnum][itemnum]
                leftitem = 0
                rightitem = 0
                if(colnum == 0):
                    leftitem = maxcols-1
                else:
                    leftitem = colnum-1
                if(colnum == maxcols-1):
                    rightitem = 0
                else:
                    rightitem=colnum+1
                leftoverhang = ""
                rightoverhang = ""
                leftvalue = self.outitems[leftitem][itemnum].value
                rightvalue = self.outitems[rightitem][itemnum].value
                logiclist = np.array([True]*len(self.p))
                if(leftvalue!=""):
                    try:
                        leftoverhang=self.p[self.p.part == leftvalue].right.iloc[0]
                    except IndexError:
                        #this means we didn't find the part!
                        raise ValueError("part {} has incorrect right overhang!".format(leftvalue))
                    if((self.outitems[-3][itemnum].value!='gibson') \
                                                and ('UNS' in leftoverhang)):
                        pass
                    else:
                        logiclist &= (self.p.left==leftoverhang)

                    #print(leftoverhang)
                if(rightvalue!=""):
                    try:
                        rightoverhang=self.p[self.p.part == rightvalue].left.iloc[0]
                    except IndexError:
                        raise ValueError("part {} has incorrect right overhang!".format(rightvalue))
                    if((self.outitems[-3][itemnum].value!='gibson') \
                                                and ('UNS' in rightoverhang)):
                        pass
                    else:
                        logiclist &= (self.p.right==rightoverhang)
                    #print(rightoverhang)
                #print("this part wants {} and {}".format(leftoverhang,rightoverhang))
                self.holdup=True
                prevval = curitem.value
                oplist,newval = self.generateOptionsList(self.p[logiclist],\
                                self.outitems[colnum][0].value,\
                                prevval,self.listEverything.value)
                curitem.options = oplist
                curitem.value = newval
        self.holdup=False

    def incrementWellPos(self,position):
        """increments a 384 well plate location such as A1 to the next logical
        position, going left to right, top to bottom"""
        poslet = self.PlateLetters.index(position[0])
        posnum = int(position[1:])
        newposlet = poslet
        newposnum = posnum+1
        if(newposnum > self.PlateRowsCols[1]):
            newposnum-=self.PlateRowsCols[1]
            newposlet+=1
        newposition = self.PlateLetters[newposlet]+str(newposnum)
        return newposition
    def finishAndSave(self,b):
        outfiletext = ",".join(self.collabels[:-1]+["targwell"])+"\n"
        outfname = self.fname1.value+".csv"
        startPos = self.DestWell.value
        curpos = startPos
        for i in range(len(self.outitems[0]))[1:]:
            outlst = []
            for nam,col in zip(self.collabels,self.outitems):
                if(nam != ""):
                    outlst+=[col[i].value]
            outlst+=[curpos]
            curpos = self.incrementWellPos(curpos)
            outfiletext+=",".join(outlst)+"\n"
        with open(os.path.join(self.mypath,"assemblies",outfname),"w") as outfle:
            outfle.write(outfiletext)
        assemfpath = os.path.join(self.mypath,"assemblies",outfname)
        #print("wrote {}".format())
        display(FileLink(assemfpath))
        display(pd.read_csv(os.path.join(self.mypath,"assemblies",outfname)))
        #b.disabled=True


    def addWidgetRow(self,labonly=True,copyrow=None):

        outcolnum=0
        for col in self.collabels:
            if(labonly):
                interwidg = widgets.Label(col)
            else:

                if(col=="name"):
                    newname = ""
                    #print(copyrow)
                    if(type(copyrow)==list):
                        newname = copyrow[outcolnum]
                    elif(type(copyrow)==int):
                        oldname = self.outitems[outcolnum][copyrow].value
                        newname = incrementString(oldname)
                    interwidg = widgets.Text(\
                            layout=self.ddlay,\
                            value=str(newname))
                elif(col==""):
                    but1 = widgets.Button(\
                        description='+',
                        button_style='success',
                        tooltip='row '+str(len(self.outitems[0])-1),
                        layout=self.eblay
                    )

                    but2 = widgets.Button(\
                        description='-',
                        button_style='danger',
                        tooltip='row '+str(len(self.outitems[0])-1),
                        layout=self.eblay,
                        #disabled=disbut
                    )
                    but1.on_click(self.add_row)
                    but2.on_click(self.remove_row)
                    interwidg =widgets.HBox([but1,but2])
                else:
                    oplist = []
                    prevval = ""
                    if(type(copyrow)==int):
                        prevval = self.outitems[outcolnum][copyrow].value
                    elif(type(copyrow)==list):
                        prevval = copyrow[outcolnum]
                    oplist, prevval = self.generateOptionsList(self.p,col,\
                                            prevval,self.listEverything.value)
                    #print(oplist)
                    #print("value is")
                    #print(prevval)
                    interwidg = widgets.Dropdown(\
                            options=oplist,\
                            value=prevval,\
                            layout=self.ddlay)
                    interwidg.observe(self.updatePartOptions,names='value')

            try:
                self.outitems[outcolnum]+=[interwidg]
            except IndexError:
                self.outitems+=[[interwidg]]
            outcolnum +=1
        self.updatePartOptions()
        for a in self.outitems[-1]:
            try:
                if(len(self.outitems[0])<=2):
                    a.children[1].disabled=True
                else:
                    a.children[1].disabled=False
            except AttributeError:
                pass


def make_assembly_file(mypath=".",externalDF = None):
    """this function will assist the user with making assembly .csv files!"""
    x=assemblyFileMaker(mypath=mypath,partsdf=externalDF)

def process_assembly_file(mypath=".",printstuff=True,partsdf=None,annotateDF=None):
    oplist = findFilesDict(os.path.join(mypath,"assemblies"))
    if(type(partsdf)==pd.DataFrame):
        parts = {"google doc":"google doc"}
    else:
        parts = findPartsListsDict(os.path.join(mypath,"partslist"))

    drop1 = widgets.Dropdown(
        options=oplist,
        #value=2,
        description='Assembly:',
    )
    drop2 = widgets.Dropdown(
        options=parts,
        #value=2,
        description='parts list:',
    )
    but = widgets.Button(
        description='Select',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        #icon='check'
    )

    #button = widgets.Button(description="Click Me!")
    #display(button)
    #print(oplist)
    def on_button_clicked(b):
        pbar = widgets.FloatProgress(
            min=0,
            max=1.0
        )
        display(pbar)
        if(drop1.value[-4:]=="xlsx" or drop1.value[-3:]=="xls"):
            x=pd.read_excel(drop1.value)
        else:
            x=pd.read_csv(drop1.value)
        if(type(partsdf)==pd.DataFrame):
            p = partsdf
        else:
            dfs = pd.read_excel(drop2.value,None)
            #print(drop1.value)

            sheetlist = list(dfs.keys())
            p = pd.DataFrame.append(dfs["parts_1"],dfs["Gibson"])

        makeEchoFile(p,x,fname = drop1.value, \
                    output = os.path.join(mypath,"output","output.csv"),\
                    sepfilename=os.path.join(mypath,"output","outputLDV.csv"),\
                    printstuff=printstuff,progbar=pbar,mypath=mypath,annotateDF=annotateDF)

        #print(drop1.value+" and "+drop2.value)

    but.on_click(on_button_clicked)
    cbox = widgets.HBox([drop1,drop2,but])
    display(cbox)

#def fixPart(partseq,enz="BsaI",circ=True,end5p=0,end3p=0,goodends=ENDDICT):

def drawConstruct(ax,construct,dnaline=3,dnascale=2,annotateDF=None,schematic=True,labels='off',showscars=0):
    """creates a dnaplotlib image of a construct in dnaseqrecord format!"""
    def substring_indexes(substring, string):
        """
        Generate indices of where substring begins in string

        >>> list(find_substring('me', "The cat says meow, meow"))
        [13, 19]
        """
        last_found = -1  # Begin at -1 so the next position to search from is 0
        while True:
            # Find next index of substring, by starting after its last known position
            last_found = string.find(substring, last_found + 1)
            if last_found == -1:
                break  # All occurrences have been found
            yield last_found

    dr = dpl.DNARenderer(scale = dnascale,linewidth=dnaline)
    part_renderers = dr.SBOL_part_renderers()

    conlist = []
    if(type(annotateDF)==pd.DataFrame):
        str_conseq = str(construct.seq).lower()
        #print("annotating!")
        #now we annotate the plasmid!!
        for feature_index in annotateDF.index:
            fname = annotateDF.iloc[feature_index]["name"]
            #iterate through all the features and see if they are in our sequence
            #but the problem is that it could be circular
            featseq = annotateDF.iloc[feature_index].sequence.lower()

            colorstr = annotateDF.iloc[feature_index].colorlist
            colorstr2 = annotateDF.iloc[feature_index].colorlist2

            #print(featcolor)
            feattype = annotateDF.iloc[feature_index].type
            featlen = len(featseq)
            #print(featcolor)
            if(featseq[-3:]=="..."):
                featseq=featseq[:-3]
            rcfeatseq = str(Dseq(featseq).rc()).lower()
            #if(feattype == 'CDS'):
                #print(featseq[:10]+"..."+featseq[-10:])
            if(featseq in str_conseq):
                #it could be in there multiple times

                for featfound in substring_indexes(featseq,str_conseq):
                    #every time we find the feature...
                    construct.add_feature(featfound,featfound+featlen,seq=None,type=feattype,label=fname,strand=1 )
                    construct.features[-1].qualifiers["color"]=colorstr
                    construct.features[-1].qualifiers["color2"]=colorstr2
            if(rcfeatseq in str_conseq):
                for featfound in substring_indexes(rcfeatseq,str_conseq):
                    #every time we find the feature...
                    construct.add_feature(featfound,featfound+featlen,seq=None,type=feattype,label=fname ,strand=-1)
                    construct.features[-1].qualifiers["color"]=colorstr
                    construct.features[-1].qualifiers["color2"]=colorstr2

    if(schematic==False):
        seqlen = len(construct)
        sp = {'type':'EmptySpace', 'name':'base', 'fwd':True, \
                                            'opts':{'x_extent':seqlen+10}}
        design = [sp]
        start,end = dr.renderDNA(ax,design,part_renderers)
    sbol_featlist = []
    flist = sorted(construct.features,key=lambda a: a.location.start)
    for feature in flist:
        #feature = a[1]
        featname = feature.qualifiers["label"]
        feattype = feature.type
        if("color" in feature.qualifiers):
            colorstr = feature.qualifiers["color"]
            if(colorstr != "(255,255,255)" and not type(colorstr)==float):
                #don't add pure white as a color
                featcolor = tuple([float(a)/255.0 for a in colorstr[1:-1].split(",")])
            else:
                featcolor = None
        else:
            colorstr = None
            featcolor = None
        if("color2" in feature.qualifiers):

            colorstr2 = feature.qualifiers["color2"]
            if(colorstr2 != "(255,255,255)" and not type(colorstr2)==float):
                #don't add pure white as a color
                featcolor2 = tuple([float(a)/255.0 for a in colorstr2[1:-1].split(",")])
            else:
                featcolor2 = None
        else:
            colorstr2 = None
            featcolor2 = None

        #print(featcolor)
        #print(feature.location)
        loclist = [feature.location.start,feature.location.end]
        if(loclist[1]<loclist[0]):
            featstrand = False
        else:
            featstrand = True
        if(feature.strand==-1):
            featstrand = False
        featstart = min(loclist)
        featend = max(loclist)
        featlen = featend-featstart
        if(not schematic):
            feat = {'type':feattype, 'name':featname, 'fwd':featstrand, \
                                    'start':featstart,'end':featend,\
                                    'opts':{'label':featname,'label_size':13,\
                                    'label_y_offset':-5,'x_extent':featlen}}
        else:
            feat = {'type':feattype, 'name':featname, 'fwd':featstrand, \
                                    #'start':featstart,'end':featend,\
                                    'opts':{'label':featname,'label_size':13,\
                                    'label_y_offset':-5}}
            if(feattype == 'CDS'):
                feat['opts']['x_extent']=30
            if(not (featcolor == None) ):
                #only add the color if it exists
                feat['opts']['color']=featcolor
            if(not (featcolor2 == None) ):
                #only add the color if it exists
                feat['opts']['color2']=featcolor2
        if(labels=="off"):
            feat['opts']['label']=""
        if(feattype == 'Scar' and not showscars):
            pass
        else:
            sbol_featlist+=[feat]

    if(schematic):
        start,end = dr.renderDNA(ax,sbol_featlist,part_renderers)
    else:
        for feat in sbol_featlist:
            dr.annotate(ax,part_renderers,feat)
    if(not construct.linear):
        vheight = 5
        curves = (end-start)*.05
        plasmid = FancyBboxPatch((start-curves, -vheight*2), \
                            (end-start)+(end-start)*.1+curves*2, vheight*2,\
                fc="none",ec="black", linewidth=dnaline, \
                boxstyle='round,pad=0,rounding_size={}'.format(curves), \
                joinstyle="round", capstyle='round',mutation_aspect=vheight/curves)
        ax.add_patch(plasmid)
    else:
        curves = 0
    ax.set_xlim([start-1.2*curves, end+1.2*curves+(end-start)*.1*(1-construct.linear)])
    ax.set_ylim([-12,12])
    #ax_dna.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
def runProgram():
    """runs the process_assembly_file function with command line prompts.
    Probably doesn't work"""
    #x=pd.read_csv(insheet,sep=",")
    #pickhand = raw_input("is this for the echo? (y/n)")
    pickhand = 'y'
    xl_file=pickPartsList()
    x,fname=pickAssembly()
    #enz=pickEnzyme()
    #p=pd.read_csv("partslist/CIDAR_parts_plate_ASS.csv",sep=",")


    #pd.ExcelFile("partslist/CIDAR_parts_plate_ASS.xlsx")
    dfs = {sheet_name: xl_file.parse(sheet_name)
          for sheet_name in xl_file.sheet_names}
    sheetlist = list(dfs.keys())
    p = pd.DataFrame.append(dfs["parts_1"],dfs["Gibson"])
    #print(p)
    try:
        if(pickhand=="n"):
            makeHandFile(p,x)
        else:
            makeEchoFile(p,x,fname = drop1.value, \
                        output = os.path.join(".","output","output.csv"),\
                        sepfilename=os.path.join(".","output","outputLDV.csv"))
    except ValueError as error:
        print("=========ERROR========")
        print(error)
        print ("")
        return False
    return True

if(__name__=="__main__"):
    progran = False
    while(not progran):
        progran = runProgram()
        if( not progran):
            z=input("Press Enter to try again")
        print("======================")
        print("")
        print("")
    z=input("Done! Press Enter to exit")




#def askForParts():
#    print "Searching for compatible parts lists..."
