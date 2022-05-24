import MDAnalysis
from scipy.spatial import distance_matrix
import sys

# CUTOFF for SB
CUT_=3.2

# load universe
u = MDAnalysis.Universe(sys.argv[1], sys.argv[2])
# number of frames
nframes = int(len(u.trajectory))

# salt bridges dictionary
# RCOO- group of ASP (D): OD1, OD2 or GLU (E): OE1 OE2
# RHN3+ from LYS (K):NZ or (RNHC(NH2)2+) of ARG (R): NH1, NH2

# list of dictionaries of negatively charged residues, divided by chain
negres=[]
# list of dictionaries of positively charged residues, divided by chain
posres=[]

# cycle on chains
for ch in ["A","B"]:
    # negatively charged residues
    dict_tmp = {}
    # select ASP atoms 
    sel = "resname ASP and name OD1 OD2 and segid "+ch 
    at = u.select_atoms(sel)
    # divide by residue 
    for r in at.residues:
        dict_tmp[("D",r.resid,ch)] = at.select_atoms("resid "+str(r.resid))
    # select GLU atoms 
    sel = "resname GLU and name OE1 OE2 and segid "+ch
    at = u.select_atoms(sel)
    # divide by residue 
    for r in at.residues:
        dict_tmp[("E",r.resid,ch)] = at.select_atoms("resid "+str(r.resid))
    # add dictionary to list
    negres.append(dict_tmp)
    # positively charged residues
    dict_tmp = {}
    # select LYS atoms
    sel = "resname LYS and name NZ and segid "+ch 
    at = u.select_atoms(sel)
    # divide by residue 
    for r in at.residues:
        dict_tmp[("K",r.resid,ch)] = at.select_atoms("resid "+str(r.resid))
    # select ARG atoms
    sel = "resname ARG and name NH1 NH2 and segid "+ch 
    at = u.select_atoms(sel)
    # divide by residue 
    for r in at.residues:
        dict_tmp[("R",r.resid,ch)] = at.select_atoms("resid "+str(r.resid))
    # add dictionary to list
    posres.append(dict_tmp)

# prepare SB dict
SB_dict={}

# cycle over frames
for i in range(0, nframes):
    # load frame
    u.trajectory[i]
    # cycle on possible pairs
    for key1 in negres[0]:
     for key2 in posres[1]:
         mdist = distance_matrix(negres[0][key1].positions,posres[1][key2].positions).min()
         # if salt-bridge is formed, add to dictionary
         if(mdist<=CUT_):
            key = (key1,key2)
            if(key in SB_dict): SB_dict[key] += 1.0
            else:               SB_dict[key]  = 1.0
    # cycle on possible pairs
    for key1 in posres[0]:
     for key2 in negres[1]:
         mdist = distance_matrix(posres[0][key1].positions,negres[1][key2].positions).min()
         # if salt-bridge is formed, add to dictionary
         if(mdist<=CUT_):
            key = (key1,key2)
            if(key in SB_dict): SB_dict[key] += 1.0
            else:               SB_dict[key]  = 1.0

# print final table to file
log=open("sbridge.dat","w")
for key in SB_dict:
    log.write("%1s %d %1s * %1s %d %1s  %lf\n" % (key[0][0],key[0][1],key[0][2],key[1][0],key[1][1],key[1][2],SB_dict[key]/float(nframes)))
