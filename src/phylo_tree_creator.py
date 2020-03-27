
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib


class PhyloTreeCreator:
    aligned_file = None
    searcher = None
    descriptions = {}
    distance_matrix = None

    def __init__(self, input_filename, searcher):
        # Read the aligned sequences and align
        self.aligned_file = AlignIO.read(input_filename, format='clustal')
        self.searcher = searcher

        # Calculate the distance matrix
        calculator = DistanceCalculator('identity')
        self.distance_matrix = calculator.get_distance(self.aligned_file)

        matplotlib.rc('font', size=6)

    def display(self):
        # Create description to be shown on the tree
        self.create_description_labels()

        # Print the distance Matrix
        print('\nDistance Matrix\n===================')
        print(self.distance_matrix)

        # Construct the phylogenetic tree using UPGMA algorithm
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(self.distance_matrix)

        self.draw_tree(tree)

    def create_description_labels(self):
        for record in self.aligned_file:
            sequence = self.searcher.get_sequence_from_genbank(record.id)
            try:
                self.descriptions[sequence.id] = self.searcher.extract_metadata(sequence)
                print(self.descriptions[sequence.id])
            except KeyError as e:
                import traceback
                traceback.print_exc()
                pass


class PhyloTreePerSpecie(PhyloTreeCreator):

    def draw_tree(self, tree):
        # Draw the phylogenetic tree
        Phylo.draw(tree, label_func=self.label_callback_species)

    def label_callback_species(self, node):
        try:
            metadata = self.descriptions[node.name]
            if 'host' in metadata:
                return str(metadata['host'])
            elif 'organism' in metadata:
                return str(metadata['organism'])
            else:
                return str(metadata['isolate_or_strain'])
        except KeyError:
            return ''

class PhyloTreePerCountry(PhyloTreeCreator):

    def draw_tree(self, tree):
        # Draw the phylogenetic tree
        Phylo.draw(tree, label_func=self.label_callback_per_country)

    def label_callback_per_country(self, node):
        try:
            metadata = self.descriptions[node.name]
            return metadata['country']
        except KeyError:
            return ''