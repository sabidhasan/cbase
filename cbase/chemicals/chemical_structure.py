import re
import pysmiles
import networkx as nx
from molmass import Formula

# Constants
ELECTRONEGATIVE_ATOMS = ('N', 'O', 'F')
HYDROGEN = 'H'

class ChemicalStructure():
  """
    Represents a chemical structure of a SMILES formula
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
      substructure_smiles             SMILES string for substructure to compare this against
    """

    substructure_graph = pysmiles.read_smiles(substructure_smiles)

    # Determine all possible candidates that map the substructure to this graph's topology
    graph_matcher = nx.isomorphism.GraphMatcher(self.graph, substructure_graph)
    candidate_mappings = graph_matcher.subgraph_isomorphisms_iter()

    # The candidates ignore atom type (nodes) and bond types (edges), so do that manually
    for candidate_dict in candidate_mappings:
      # Each dict maps nodes from self.graph to `other` graph
      if self._candidate_is_valid(candidate_dict, substructure_graph):
        return True

    return False


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


  def _candidate_is_valid(self, candidate_dict, substructure_graph):
    """
        Determines if the candidate represented by candidate_dict is valid for the desired
        substructure - this means atom types and bond types match
        
        candidate_dict            dict - maps nodes in self.graph to those in a substructure
        substructure_graph        networkx object
    """
    
    # Candidate dict maps node IDs in self.graph to that in the substructure - this the reverse
    substructure_id_to_parent_dict = {v: k for k, v in candidate_dict.items()}

    # A candidate mapping is valid if all atoms have same type, and all bond orders are valid
    for parent_atom_id, candidate_atom_id in candidate_dict.items():
      # Get the atom descriptions out of the graph from the atom IDs
      parent_atom = self.graph.nodes[parent_atom_id]
      candidate_atom = substructure_graph.nodes[candidate_atom_id]
      if not self._atoms_are_identical(parent_atom, candidate_atom):
          return False

      # Check that all bonds from this atom match the parent's
      bonds_from_candidate_atom = substructure_graph.adj[candidate_atom_id].items()
      if not self._bonds_are_identical(bonds_from_candidate_atom, substructure_id_to_parent_dict, parent_atom_id):
        return False

    return True


  @staticmethod
  def _atoms_are_identical(element_dict_1, element_dict_2):
    # If the atom types are different, this dict is not valid
    return (element_dict_1['element'] == element_dict_2['element'] and element_dict_1['charge'] 
            == element_dict_2['charge'] and element_dict_1['aromatic'] == element_dict_2['aromatic'])


  def _bonds_are_identical(self, bonds_from_candidate_atom, substructure_id_to_parent_dict, parent_atom_id):
    # bonds_from_substructure_atom 
    for substructure_bond_destination, substructure_bond_properties in bonds_from_candidate_atom:
      # get bond order for bond from parent_atom_id to equivalent of `substructure_bond_destination`
      parent_destination_bond_id = substructure_id_to_parent_dict[substructure_bond_destination]
      parent_bond_order = self.graph.get_edge_data(parent_atom_id, parent_destination_bond_id)['order']
      if substructure_bond_properties['order'] > parent_bond_order:
        return False

    return True
