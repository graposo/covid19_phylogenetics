
from Bio import Entrez
from Bio import SeqIO
from datetime import datetime

class SearcherConfiguration:
    email = "A.N.Other@example.com"
    species_max = 10
    species_sequences_filename = None
    countries_max = 15
    countries_sequences_filename = None
    total_search_max = 500


class CoronaVirusSearcher:
    configuration = None
    species = {}
    countries = {}
    date_cov2_baseline = datetime.strptime('01-Dec-2019', '%d-%b-%Y')

    def __init__(self, configuration):
        self.configuration = configuration
        Entrez.email = configuration.email

        if self.configuration.species_sequences_filename is None:
            self.configuration.species_sequences_filename = 'species_sequences_%d.fasta' % self.configuration.total_search_max
            configuration.species_sequences_filename = self.configuration.species_sequences_filename

        if self.configuration.countries_sequences_filename is None:
            self.configuration.countries_sequences_filename = 'countries_sequences_%d.fasta' % self.configuration.total_search_max
            configuration.countries_sequences_filename = self.configuration.countries_sequences_filename


    def find_and_filter(self):
        # Download sequences from DB and filter them
        self.__find_sequences(self.__include_in_result, '(coronavirus[organism]) AND genome[title]')

        # Write sequences to disk
        self.__write_sequences_to_file(self.species.values(), self.configuration.species_sequences_filename)
        self.__write_sequences_to_file(self.countries.values(), self.configuration.countries_sequences_filename)

        return self.species, self.countries

    def extract_metadata(self, sequence):
        result = {'description': sequence.annotations['source'],
                  'date': sequence.annotations['date'],
                  'country': sequence.features[0].qualifiers['country'][0],
                  'collection_date': sequence.features[0].qualifiers['collection_date'][0],
                  'organism': sequence.features[0].qualifiers['organism'][0]}

        if 'host' in sequence.features[0].qualifiers:
            result['host'] = sequence.features[0].qualifiers['host'][0]
        else:
            # For non Homo sapiens there seems to not be host metadata, but the organism field can give us a hint of of what animal is it
            result['host'] = result['organism']

        if 'isolate' in sequence.features[0].qualifiers:
            result['isolate_or_strain'] = sequence.features[0].qualifiers['isolate'][0]
        elif 'strain' in sequence.features[0].qualifiers:
            result['isolate_or_strain'] = sequence.features[0].qualifiers['strain'][0]
        else:
            raise KeyError('Not found isolate neither strain fields')

        return result

    def __is_bat_sequence(self, sequence_metadata):
        return 'host' in sequence_metadata and 'bat' in sequence_metadata['host'].lower()

    def __is_human_sequence(self, sequence_metadata):
        return 'host' in sequence_metadata and 'sapiens' in sequence_metadata['host'].lower()


    def __write_sequences_to_file(self, sequences, filename):
        with open(filename, 'w') as output_file:
            SeqIO.write(sequences, filename, format='fasta')

    def __find_sequences(self, filtering_callback, search_string):

        result = []
        with Entrez.esearch(db="nucleotide", retmax=self.configuration.total_search_max,
                            term=search_string) as search_handle:
            # Iterate over the ids and retieve the genbank record
            search_record = Entrez.read(search_handle)
            for seq_id in search_record['IdList']:
                sequence_record = self.get_sequence_from_genbank(seq_id)
                if filtering_callback(sequence_record):
                    return result

        return result


    def __include_in_result(self, sequence_record):
        try:
            sequence_metadata = self.extract_metadata(sequence_record)
        except KeyError:
            return

        # Make sure to have only one by specie
        if not self.__is_human_sequence(sequence_metadata) and sequence_metadata['host'] not in self.species.keys() and len(self.species.keys()) < self.configuration.species_max:
            print('Adding species [%s]' % sequence_metadata['host'])
            self.species[sequence_metadata['host']] = sequence_record
        elif self.__is_human_sequence(sequence_metadata) and not 'homo_sapiens' in self.species.keys():     # Always add one homo sapiens
            print('Adding species [%s]' % 'homo_sapiens')
            self.species['homo_sapiens'] = sequence_record

        # Make sure to have only 1 per country
        try:
            collection_date = sequence_metadata['collection_date']
            if self.__is_human_sequence(sequence_metadata) and self.__after_cov_baseline_time(collection_date) and \
                    sequence_metadata['country'] not in self.countries.keys() and len(self.countries.keys()) < self.configuration.countries_max:
                print('Adding country [%s] with collection_date [%s] for host [%s]' % (sequence_metadata['country'], collection_date, sequence_metadata['host']))
                self.countries[sequence_metadata['country']] = sequence_record
            elif self.__is_bat_sequence(sequence_metadata) and not 'bat' in self.countries.keys() and \
                     len(self.countries.keys()) < self.configuration.countries_max:
                print('Adding species [%s] to the countries dataset' % sequence_metadata['host'])
                self.countries['bat'] = sequence_record
        except ValueError:
            # Exception will be thrown if 'collection_date' field has a different date format. We just ignore the record
            pass

        # Return true if there is no need to keep downloading sequences
        return  len(self.species.keys()) >= self.configuration.species_max and \
                len(self.countries.keys()) >= self.configuration.countries_max


    def __after_cov_baseline_time(self, collection_date):
        return datetime.strptime(collection_date, '%d-%b-%Y') > self.date_cov2_baseline


    def get_sequence_from_genbank(self, seq, debug=False):
        with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=seq) as handle:
            seq_record = SeqIO.read(handle, 'gb')

        return seq_record


