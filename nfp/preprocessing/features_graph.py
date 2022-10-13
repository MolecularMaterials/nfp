import networkx as nx
from mendeleev import element

class Tokenizer(object):
    """ A class to turn arbitrary inputs into integer classes. """
    
    def __init__(self):
        # the default class for an unseen entry during test-time
        self._data = {'unk': 1}
        self.num_classes = 1
        self.train = True
        self.unknown = []
    
    def __call__(self, item):
        """ Check to see if the Tokenizer has seen `item` before, and if so,
        return the integer class associated with it. Otherwise, if we're
        training, create a new integer class, otherwise return the 'unknown'
        class.

        """
        try:
            return self._data[item]

        except KeyError:
            if self.train:
                self._add_token(item)
                return self(item)

            else:
                # Record the unknown item, then return the unknown label
                self.unknown += [item]
                return self._data['unk']
                
    def _add_token(self, item):
        self.num_classes += 1
        self._data[item] = self.num_classes

# The rest of the methods in this module are specific functions for computing
# atom and bond features. New ones can be easily added though, and these are
# passed directly to the Preprocessor class.

def atom_features_v1(node):
    """ Return an integer hash representing the atom type
    """
    return str((
        node['element'],
        node['coordinationNumber'],
        node['type'],
        node['neighborAtoms']
    ))

def atom_features_v2(node):
    """ Return an integer hash representing the atom type
    """
    return str((
        node['element'],
        node['coordinationNumber'],
        node['type']
    ))


def atom_features_CGCNN(node):
    return str((
    #    element(node['element']).atomic_number,
    #    node['coordinationNumber'],
    #    node['type'],
        element(node['element']).electronegativity('sanderson'),
        element(node['element']).covalent_radius_cordero,
        element(node['element']).electron_affinity,
        element(node['element']).ionenergies[1],
        element(node['element']).nvalence(),
        element(node['element']).group_id,
        element(node['element']).period,
        element(node['element']).block,
        element(node['element']).atomic_volume
    ))

    

def bond_features_v1(edge, **kwargs):
    """ Return an integer hash representing the bond type.
    
    flipped : bool
        Only valid for 'v3' version, whether to swap the begin and end atom types

    """
  #  beginAtom = ''.join([i for i in edge[0] if not i.isdigit()])
  #  endAtom = ''.join([i for i in edge[1] if not i.isdigit()])
    return str((
        edge['bondLength'],
        edge['bondType'],
        sorted([
            edge['bondAtom1'],
            edge['bondAtom2']]),
        ))

def bond_features_noBondLength(edge, **kwargs):
    return str((
        edge['bondType'],
        sorted([
            edge['bondAtom1'],
            edge['bondAtom2']]),
        ))

def bond_features_v3(bond, **kwargs):

    return str((
        bond.GetBondType(),
        bond.GetIsConjugated(),
        bond.GetStereo(),
        get_ring_size(bond),
        sorted([
            bond.GetBeginAtom().GetSymbol(),
            bond.GetEndAtom().GetSymbol()]),
        ))


def bond_features_v3(bond, flipped=False):

    if not flipped:
        start_atom = atom_features_v1(bond.GetBeginAtom())
        end_atom = atom_features_v1(bond.GetEndAtom())

    else:
        start_atom = atom_features(bond.GetEndAtom())
        end_atom = atom_features(bond.GetBeginAtom())

    return str((
        bond.GetBondType(),
        bond.GetIsConjugated(),
        bond.GetStereo(),
        get_ring_size(bond),
        bond.GetEndAtom().GetSymbol(),
        start_atom,
        end_atom))


#if __name__ == "__main__":
#    G = nx.read_graphml('test.graphml')
#    for m,edge in enumerate(G.edges):
#        if 'Mo13' in edge:
#            print(edge)
    