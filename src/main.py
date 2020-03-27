from Bio import Entrez
from Bio import SeqIO

from Bio.Align.Applications import ClustalwCommandline
from search_engine import SearcherConfiguration, CoronaVirusSearcher
from phylo_tree_creator import PhyloTreePerCountry, PhyloTreePerSpecie

def align_sequences(filename_input, filename_output):
    ccli = ClustalwCommandline('clustalw2', infile=filename_input, outfile=filename_output)
    print('Aligning sequences from file [%s]...' % filename_input)
    output, error = ccli()
    print(output)


if __name__ == '__main__':
    config = SearcherConfiguration()
    searcher = CoronaVirusSearcher(config)

    # Download sequences and store them on disk
    #searcher.find_and_filter()


    #align_sequences('src/species_sequences_500.fasta', 'species.alg')
    #tree_creator = PhyloTreePerSpecie('species.alg', searcher)

    #align_sequences('src/countries_sequences_500.fasta', 'countries.alg')
    tree_creator = PhyloTreePerCountry('countries.alg', searcher)
    tree_creator.display()