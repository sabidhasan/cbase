import re
import pysmiles
import networkx as nx
from molmass import Formula

# Constants
ELECTRONEGATIVE_ATOMS = ('N', 'O', 'F')
HYDROGEN = 'H'

class ChemicalStructure():
  """
    Represents a chemical structure
  """

  def __init__(self, smiles_string):
      """
        smiles_string           str of smiles representation of the molecule
      """

      self.smiles = smiles_string
      try:
        # NetworkX graph representation of the molecule
        self.graph = pysmiles.read_smiles(smiles_string, explicit_hydrogen=True)
      except:
        raise ValueError('Invalid structure ' + smiles_string)
  

  def __contains__(self, substructure_smiles):
      """
        Returns whether substructure occurs in this structure
      """

      def atoms_match(atom1, atom2):
        ret = (atom1['element'] == atom2['element'] and atom1['aromatic'] == atom2['aromatic'])
        print(atom1, atom2, ret)
        return ret
      
      # def bonds_match(bond1, bond2):
      #   return bond1['order'] >= bond1['order']

      substructure = pysmiles.read_smiles(substructure_smiles)
      nx.isomorphism.is_isomorphic(self.graph, substructure, node_match=atoms_match)
                                         #edge_match=bonds_match)


  @property
  def ring_count(self):
    """
      Returns number of rings in molecule
    """

    return len(nx.cycle_basis(self.graph))


  @property
  def mol_weight(self):
    """
      Returns molecular weight of molecule
    """
    
    atoms_list = nx.get_node_attributes(self.graph, 'element').values()
    return sum(Formula(atom).mass for atom in atoms_list)   

    
  @property
  def hill_formula(self):
    """
      Return Hill notaton formula for compound
    """

    atoms_list = nx.get_node_attributes(self.graph, 'element').values()
    return Formula(''.join(atoms_list)).formula


  @property
  def num_h_bonds(self):
    h_donors, h_acceptors = 0, 0

    for atom_index in self.graph.nodes:
      atom_element = self.graph.nodes[atom_index]['element']
      if atom_element == HYDROGEN:
        # Check if neighbor is N, O, F - we pick index 0, because H has valence of 1
        neighbor_index = list(self.graph.edges(atom_index))[0][1]
        neighbor_element = self.graph.nodes[neighbor_index]['element']
        
        if neighbor_element in ELECTRONEGATIVE_ATOMS:
          h_donors += 1
      elif atom_element in ELECTRONEGATIVE_ATOMS:
        h_acceptors += 1
    
    return h_donors, h_acceptors


  @property
  def heavy_atom_count(self):
    """
      Return count of heavy atoms
    """
    
    atoms_list = nx.get_node_attributes(self.graph, 'element').values()
    return len([atom for atom in atoms_list if atom != HYDROGEN])


  @staticmethod
  def _cas_number_is_valid(cas_string):
    """
      Perform basic validation for CAS number
    """

    try:
  	  return re.search(r'\d{2,7}-\d{2}-\d', cas_number_input, re.I) is not None
    except:
      return False
