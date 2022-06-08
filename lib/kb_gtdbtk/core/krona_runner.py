'''
Run KronaTools
'''

import json
from collections import Counter

# create a Krona chart based on a text file
# https://github.com/marbl/Krona/wiki/Importing-text-and-XML-data
def run_krona_import_text(
        runner: Callable[[List[str]], None],
        output_dir: Path,
        temp_dir: Path,
        ) -> None:

    input_file = _create_input_file(output_dir, temp_dir)

    chart_file_name = 'krona_chart.html' # if updated, the index.html file should also be updated accordingly
    output_file = os.path.join(output_dir, chart_file_name)

    import_text_cmd = ['ktImportText', input_file, '-o', output_file]

    logging.info('Starting Command:\n' + ' '.join(import_text_cmd))
    runner(import_text_cmd)

def _create_input_file(output_dir: Path, temp_dir: Path) -> Path:
    '''
    create tab-delimited text file from gtdbtk result files ['gtdbtk.ar53.summary.tsv.json', 'gtdbtk.bac120.summary.tsv.json']
    '''
    summary_files = ['gtdbtk.ar53.summary.tsv.json', 'gtdbtk.bac120.summary.tsv.json']

    taxonomies = list()
    for summary_file in summary_files:

        summary_file_path = os.path.join(output_dir, summary_file)

        if os.path.exists(summary_file_path):

            with open(summary_file_path, 'r') as myfile:
                data=myfile.read()

            obj = json.loads(data)

            for genome_data in obj.get('data', []):
                classification = genome_data.get('classification', '')
                taxonomies.append(classification)

    taxo_counter = Counter(taxonomies)
    krona_input_file = 'krona_input.txt'
    krona_input_path = os.path.join(temp_dir, krona_input_file)

    with open(krona_input_path, 'w') as krona_file:
        for taxonomy, count in taxo_counter.items():
            line = '{}\t{}\n'.format(count, taxonomy.replace(";", '\t'))
            print(line)
            krona_file.write(line)

    return krona_input_path
