from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from shutil import copyfile
import os

# copied from https://github.com/kbaseapps/kb_cmash/blob/master/lib/kb_cmash/utils/misc_utils.py


def load_fastas(callback_url, scratch, upa):
    '''
    '''
    dfu = DataFileUtil(callback_url)
    au = AssemblyUtil(callback_url)
    obj_data = dfu.get_objects({"object_refs": [upa]})['data'][0]
    obj_type = obj_data['info'][2]

    if 'KBaseSets.GenomeSet' in obj_type:
        upas = [gsi['ref'] for gsi in obj_data['data']['items']]
    elif 'KBaseSearch.GenomeSet' in obj_type:
        upas = [gse['ref'] for gse in obj_data['data']['elements'].values()]
    elif "KBaseGenomes.Genome" in obj_type:
        upas = [upa]
    elif "KBaseGenomes.ContigSet" in obj_type or "KBaseGenomeAnnotations.Assembly" in obj_type:
        # in this case we use the assembly file util to get the fasta file
        # file_output = os.path.join(scratch, "input_fasta.fa")
        faf = au.get_assembly_as_fasta({"ref": upa})
        return [(faf['path'], upa)]
    elif "KBaseSets.AssemblySet" in obj_type:
        fasta_paths = []
        for item_upa in obj_data['data']['items']:
            faf = au.get_assembly_as_fasta({"ref": item_upa['ref']})
            fasta_paths.append((faf['path'], item_upa['ref']))
        return fasta_paths

    data_objs = dfu.get_objects({"object_refs": upas})['data']
    assembly_upas = [d['data']['assembly_ref'] for d in data_objs]
    fasta_paths = []

    for upa in assembly_upas:
        # first we want to ge the associated assemblies.
        faf = au.get_assembly_as_fasta({"ref": upa})
        # faf = gfu.genome_features_to_fasta({"genome_ref": upa})
        fasta_paths.append((faf['path'], upa))
    return fasta_paths

def create_html_report(callback_url, scratch, workspace_name):
    '''
    '''
    output_dir = os.path.join(scratch, 'output')
    dfu = DataFileUtil(callback_url)
    report = KBaseReport(callback_url)
    copyfile(os.path.join(os.path.dirname(__file__), 'index.html'), 
             os.path.join(output_dir, 'index.html'))

    report_shock_id = dfu.file_to_shock({'file_path': output_dir,
                                        'pack': 'zip'})['shock_id']

    html_file = {
        'shock_id': report_shock_id,
        'name': 'index.html',
        'label': 'index.html',
        'description': 'HTML report for GTDBTk'
        }
    
    report_info = report.create_extended_report({
                    'direct_html_link_index': 0,
                    'html_links': [html_file],
                    'report_object_name': report_name,
                    'workspace_name': workspace_name
                })
    return {
        'report_name': report_info['name'],
        'report_ref': report_info['ref']
    }