import random
from functools import *
from operator import mul
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import math
import os

edges = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,\
             24,47,48,71,72,95,96,119,120,143,144,167,168,191,192,215,216,\
             239,240,263,264,287,288,311,312,335,336,359,\
             360,361,362,363,364,365,366,367,368,369,370,371,\
             372,373,374,375,376,377,378,379,380,381,382,383,384]
topbottom= [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,\
             360,361,362,363,364,365,366,367,368,369,370,371,\
             372,373,374,375,376,377,378,379,380,381,382,383]

def drawPlate(welldict={}):
    fig, ax = plt.subplots()
    # create 3x3 grid to plot the artists
    grid = np.mgrid[0:1:24j,.7:0:16j].reshape(2, -1,order="F").T

    patches = []
    uniqconstructs = set([welldict[a] for a in welldict])
    numconstructs=len(uniqconstructs)
    colrange=np.linspace(0,1,numconstructs)
    random.shuffle(colrange)
    concolor = {a:b for a,b in zip(uniqconstructs,range(numconstructs))}
    print(concolor)
    colors=plt.cm.rainbow([0]+colrange)

    condiv=0
    #if(constructs>0):
        #condiv=int(384/constructs)
    # add a rectangle
    color=colrange[0]
    #print(len(colors))
    clist=[]
    centroids = [(0,0)]*numconstructs
    pcount = [0]*numconstructs
    for a in range(384):
        i = -1
        try:
            i=concolor[welldict[a]]
            color= colors[i]
        except KeyError:
            color="lightgray"
        if(i>=0):
            centroids[i] = \
                    (centroids[i][0]+grid[a][0],\
                    centroids[i][1]+grid[a][1])
            pcount[i]+=1
        rect = mpatches.Rectangle(grid[a], .6/17, 1/28, \
                                ec="none", fc=color,alpha=.7)
        plt.gca().add_patch(rect)
        #patches.append(rect)
    #label(grid[1], "Rectangle")
    centroids = [(a[0]/b,a[1]/b) for a,b in zip(centroids,pcount)]
    for construct,centpt in zip(uniqconstructs,centroids):
        offset = (len(construct)*.0125)/2
        plt.gca().text(centpt[0]-offset,centpt[1],construct,\
                        bbox={'facecolor':'white','edgecolor':'white',\
                                            'alpha':0.9,'pad':2})
    #colors = np.linspace(0, 1, len(patches))
    collection = PatchCollection(patches,alpha=.3)
    #collection.set_array(np.array(colors))
    #ax.add_collection(collection)
    #ax.add_line(line)

    plt.axis('equal')
    #plt.xlim([0,1.2])
    #plt.ylim([0,.8])
    plt.axis('off')
    #plt.tight_layout()

    plt.show()
def allcomb(listoflists):
    """creates all paths through a list"""
    if(len(listoflists)==1):
        return([[a] for a in listoflists[0]])
    outlist = []
    for element in listoflists[0]:
        for lists in allcomb(listoflists[1:]):
            outlist+=[[element]+lists]
    return outlist

def multid_dilution(inducers,wells,construct,fname,blacklist=[],constructnames = [],inducernames=[],shuffle=False,wellorder="across",draw=True, mypath="."):
    outfile = "Source Plate Name,Source Plate Type,Source Well,    Sample ID,Sample Name,Sample Group,Sample Comment,Destination Plate Name,    Destination Well,Transfer Volume\n"
    suppfile = "Well,Construct,"
    suppfile+= ",".join(inducernames)
    suppfile+= "\n"

    eachline = "Source[1],384PP_AQ_BP,{},,,,,Destination[1],{},{}\n"
    rows = "ABCDEFGHIJKLMNOP"

    columns = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
    wellslist=range(384)
    if(wellorder=="down"):
        wellslist=list(list(np.reshape(np.reshape(wellslist,(16,24)),(1,-1),order="F"))[0])
    wellsused = []
    wellct = 0
    wlist={}
    for con in construct:
        curconname = constructnames[construct.index(con)]
        wellct = rows.index(con[0])*len(columns)+columns.index(int(con[1:])) #this starts the
        #print(wellct)
        #print(wellslist)

        wellind = wellslist.index(wellct)

        #well counter at the right row, in the leftmost column
        #if(wellorder=="down"):
        #    wellct = rows.index(con[0])*len(columns)+columns.index(int(con[1:]))
        controlwell = [[0]*len(inducers)]

        acomb = allcomb(inducers)
        if(not controlwell in acomb):
            #you should really do some un-induced controls!!
            acomb+=controlwell

        if(shuffle):
            random.shuffle(acomb)
        #print(acomb[:15])
        iterct = reduce(mul,[len(a) for a in inducers],1) #this is the number of wells per construct
        lensum = sum([len(a) for a in inducers])
        ict2 = 0
        conwells=[]
        for iteration in acomb:
            while(wellslist[wellind] in blacklist):
                wellind+=1#moves over the current well if it is in the blacklist
            wellct=wellslist[wellind]
            colct = wellct%len(columns) #deciphering the well number to column and rows
            rowct = int(wellct/len(columns))
            wellstr = rows[rowct]+str(columns[colct])
            #print(iteration)
            suppfile+= ",".join([wellstr,curconname]+[str(a) for a in iteration])+"\n"

            #each comination of volumes is a list
            #where each element corresponds to an inducer
            for ind_i in range(len(iteration)):
                #this goes through every inducer and
                #tells the echo to pipet it.

                #print(rowct)
                #print(colct)

                if(iteration[ind_i]>0):
                    #don't tell the echo to pipet zero nl!
                    indname = inducernames[ind_i]

                    usewell = wells[ind_i]
                    if(type(usewell)==list):
                        #if we have multiple sources, randomly pick one of them!
                        usewell = random.choice(usewell)
                    outfile += eachline.format(usewell,wellstr,iteration[ind_i])
            #even 0,0 will be marked as used. This is good!
            wellsused+=[wellct]
            conwells+=[wellct]
            wlist.update({wellct:curconname})

            wellind+=1#here is where we just go along left to right, top to bottom
        #wlist+=[conwells]
    if(draw):
        drawPlate(wlist)
    outfle = open(os.path.join(mypath,fname),"w")
    outfle.write(outfile)
    outfle.close()
    print("wrote "+os.path.join(mypath,fname))
    supfleout= open(os.path.join(mypath,"supp_"+fname),"w")
    supfleout.write(suppfile)
    supfleout.close()
    print("wrote "+os.path.join(mypath,"supp_"+fname))
    return wellsused
def multid_dilution_wrapper(inducers,constructs,fname,avoidedges=[],maxinducer=500,\
                wellvol=50,shuffle=False,wellorder="across",mypath="."):
    """this function contains some helpful pre-sets for doing multiple
    inducer sweeps in a 384 well plate.
    inducers:
    this is a list of lists that looks like this:
    [ ["name",[int(amount1),int(amount2)],["well",number of wells]]]

    "name"
    is the name of the inducer. for example, "ATC" or "ARA"

    amount1, amount2
    are the percent of the 'maxinducer' volume to add. for
    example, if you put '50', that means 50 percent of the max inducer volume
    (which is default 500), so the echo will add 250 nl of that inducer.

    "well" is the well where your inducers will start, in the source plate.
    You will need about 220 ul to populate an entire 384 well plate which does
    not fit into one source plate well. This is why we usually do many source
    plate wells. The first of these is denoted by "well", and the rest are placed
    vertically downwards in the plate. For example if you put "K1" here, you also
    need to put inducer into wells "L1","M1", and "N1", or more or less wells in
    a vertical column as determined by the number given by 'number of wells.'

    constructs:

    This is a list of strings that denote the names of the different constructs
    you are going to use. Each construct gets the same array of inducers.

    fname:

    this is the name of the file you want to make. Usually "date.csv". The program
    also makes "supp_<fname>" as a supplementary data file used with the biotek analysis
    package.

    avoidedges:

    a list of numbers indicating which wells to avoid. Use the two predifined lists
    called "edges" or "topbottom" to block out all edge wells or the top and bottom rows

    maxinducer:

    "100%" volume of inducer to add. Usually 500 nl so it's 100x diluted into the
    50 uL culture volume.

    wellvol:

    culture volume

    shuffle:

    True means that the wells which recieve different inducer concentrations Will
    be randomly distributed within the set of wells set aside for each construct.

    wellorder:

    "across" means lettered rows from left to right. A1 -> A2 -> A3 etc
    "down" means columns down (letters in order) then left to right. So, A1 -> B1
    ->->A2 -> B2 etc.
    """
    if(avoidedges==edges):
        #print("edges")
        start="B2"
        rows=14
        cols=22
    elif(avoidedges==topbottom):
        start="B1"
        rows=14
        cols=24
    else:
        start="A1"
        rows=16
        cols=24
    #divideplate=int(rows/len(constructs))
    conwells=[]
    if(wellorder=="across"):
        divideplate=math.ceil(rows/len(constructs))
        while(divideplate*len(constructs)>rows):
            divideplate-=1
        for i in range(len(constructs)):
            conwells+=[chr(ord(start[0])+divideplate*i)+start[1:]]
    else:
        divideplate=math.ceil(cols/len(constructs))
        while(divideplate*len(constructs)>cols):
            divideplate-=1
        for i in range(len(constructs)):
            conwells+=[start[0]+str(int(start[1:])+divideplate*i)]

    indonly=[a[0] for a in inducers]
    volonly=[[int(b*maxinducer/100) for b in a[1]] for a in inducers]
    sourceonly=[]
    for src in [a[2] for a in inducers]:
        allwells=[]
        for i in range(src[1]):
            allwells+=[chr(ord(src[0][0])+i)+src[0][1:]]
        sourceonly+=[allwells]
    print("load inducers into: ")
    for ind,amts in zip(indonly,sourceonly):
        print("{} : {}".format(ind,str(amts)))
    print()
    inducervol=50*len(sourceonly[0])*1.1
    print("make %02f ul of 100x each inducer"%(inducervol))
    print()
    print("constructs start at: ")
    for cname,cwell in zip(constructs,conwells):
        print("{} : {}".format(cname,str(cwell)))
    print("prepare {} ml of each construct for {} well volume".format(int(wellvol*1.1*384/len(constructs))/1000,wellvol))
    return multid_dilution(volonly,sourceonly,conwells,fname,avoidedges,constructs,indonly,shuffle,wellorder=wellorder,mypath=mypath)
