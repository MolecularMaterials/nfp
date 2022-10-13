from ase import neighborlist
import numpy as np
import networkx as nx
import networkx as nx
import time

start_time = time.time()

def adsGraphGen(atoms,adsorbateIndices,structName=None,totalE=None): 
    """
    Create a graph for specified adsorbates and surface atoms

    atoms(ase.Atoms object): structure
    adsorbateIndices(list): List of adsorbate indices
    """
    cutOff = np.asarray(neighborlist.natural_cutoffs(atoms))
    nl = neighborlist.NeighborList(cutOff, skin=0.25, self_interaction=False, bothways=True)
    nl.update(atoms)
    nlist = []
    coordNumList = []
    for i, atom in enumerate(atoms):
        indices, offset = nl.get_neighbors(i)
        indices = list(set(indices))
        nlist.append(indices)
        count = 0
        for i in indices:
            if i in adsorbateIndices:
                count+=1
        coordNumList.append(len(indices)-count) # bonds between adsorbed atoms are not counted as a coordination
    fsnList = [] # first coordination shell list
    fsnIndices = []
    for ads in adsorbateIndices:
        fsnList.append(nlist[ads])
        for i in nlist[ads]:
            if i not in fsnIndices and i not in adsorbateIndices:
                fsnIndices.append(i)
    ssnList = []
    for fsn in fsnList:
        ssn = [nlist[idx] for idx in fsn]
        ssnList.append(ssn)

    ssnIndices = []
    for ssn in ssnList:
        for ss in ssn:
            for i in ss:
                if i not in ssnIndices and i not in fsnIndices and i not in adsorbateIndices:
                    ssnIndices.append(i)


    tsnList = []
    for ssnList2 in ssnList:
        for ssn in ssnList2:
            tsn = [nlist[idx] for idx in ssn]
            tsnList.append(tsn)
    tsnIndices = []
    for tsn in tsnList:
        for ts in tsn:
            for i in ts:
                if i not in tsnIndices and i not in ssnIndices and i not in fsnIndices and i not in adsorbateIndices:
                    tsnIndices.append(i)


    bondList = []

    for i, ads in enumerate(adsorbateIndices): # Create adsorbate-adsorbate and adsorbate-(nearest neighbor) surface atom bonds
        for j in fsnList[i]:
            if [ads,j] not in bondList:
                bondList.append([ads,j])

    for i, fsn in enumerate(fsnList): # Create bonds between first and second shell
        for j, fs in enumerate(fsn):
            for z in ssnList[i][j]:
                if [fs,z] not in bondList:
                    bondList.append([fs,z])

    for i in ssnIndices:       # Creat bonds between atoms in the second shell
        for j in ssnIndices:
            if i!=j and i in nlist[j] and [i,j] not in bondList:
                bondList.append([i,j])


    for ssnList2 in ssnList: # Create bonds between second and third shell
        for i, ssn in enumerate(ssnList2):
            for j, ss in enumerate(ssn):
                for z in tsnList[i][j]:
                    if [ss,z] not in bondList:
                        bondList.append([ss,z])

    for i in tsnIndices:       # Creat bonds between atoms in the third shell
        for j in tsnIndices:
            if i!=j and i in nlist[j] and [i,j] not in bondList:
                bondList.append([i,j])

    finBondList = []
    dupList = []

    for i, bl1 in enumerate(bondList):
        for j, bl2 in enumerate(bondList):
            if set(bl1)==set(bl2) and i<j:
                dupList.append(bl2)

    for bl in bondList:
        if bl not in dupList:
            finBondList.append(bl)

    graph = nx.Graph(totalEnergy=totalE)
    graph.graph['structureID']=structName
    for bl in finBondList:
        a1idx = bl[0]
        a2idx = bl[1]
        atom1 = atoms[a1idx].symbol
        atom2 = atoms[a2idx].symbol
        a1lbl = atom1+str(a1idx)
        a2lbl = atom2+str(a2idx)
        d = atoms.get_distance(a1idx,a2idx,mic=True)
        graph.add_edge(a1lbl,a2lbl,bondLength=round(d,3))
        if a1idx in adsorbateIndices and a2idx in adsorbateIndices:
            graph.edges[(a1lbl,a2lbl)]['bondType']='molecular'
        elif a1idx not in adsorbateIndices and a2idx not in adsorbateIndices:
            graph.edges[(a1lbl,a2lbl)]['bondType']='lattice'
        else:
            graph.edges[(a1lbl,a2lbl)]['bondType']='adsorption'
        graph.edges[(a1lbl,a2lbl)]['bondAtom1']=atom1
        graph.edges[(a1lbl,a2lbl)]['bondAtom2']=atom2
        graph.nodes[a1lbl]['element']= atom1
        graph.nodes[a1lbl]['degree']= graph.degree[a1lbl]
        graph.nodes[a2lbl]['element']= atom2
        graph.nodes[a2lbl]['degree']= graph.degree[a2lbl]
        if a1idx in adsorbateIndices:
            graph.nodes[a1lbl]['type']='adsorbate'
        elif a1idx in fsnIndices:
            graph.nodes[a1lbl]['type']='firstCoordinationShell'
        elif a1idx in ssnIndices:
            graph.nodes[a1lbl]['type']='secondCoordinationShell'
        elif a1idx in tsnIndices:
            graph.nodes[a1lbl]['type']='thirdCoordinationShell'
        else:
            graph.nodes[a1lbl]['type']='unknown'
        graph.nodes[a1lbl]['coordinationNumber']=coordNumList[a1idx]
        graph.nodes[a1lbl]['neighborAtoms'] = str(sorted(atoms[i].symbol for i in nlist[a1idx]))

        if a2idx in adsorbateIndices:
            graph.nodes[a2lbl]['type']='adsorbate'
        elif a2idx in fsnIndices:
            graph.nodes[a2lbl]['type']='firstCoordinationShell'
        elif a2idx in ssnIndices:
            graph.nodes[a2lbl]['type']='secondCoordinationShell'
        elif a2idx in tsnIndices:
            graph.nodes[a2lbl]['type']='thirdCoordinationShell'
        else:
            graph.nodes[a2lbl]['type']='unknown'
        graph.nodes[a2lbl]['coordinationNumber']=coordNumList[a2idx]
        graph.nodes[a2lbl]['neighborAtoms'] = str(sorted(atoms[i].symbol for i in nlist[a2idx]))

    return graph