import pandas as pd
import os
import uuid
import logging
import json
import re
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenericsAPIClient import GenericsAPI

# Subsetting Class
class Subsetting_Matrices:

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.scratch = config['scratch']

        self.dfu = DataFileUtil(self.callback_url)

        # set up directory for files folder
        self.output_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(self.output_dir)
        self.files_folder = os.path.join(self.output_dir, 'files')
        os.mkdir(self.files_folder)

        self.file_paths = []
        self.html_paths = []

        self.GenAPI = GenericsAPI(self.callback_url)

    def _get_df(self, params):
        """
        Get Amplicon Matrix Data then make Pandas.DataFrame(),
        also get taxonomy data and add it to df.
        """

        logging.info('Getting DataObject')

        # Amplicon data
        obj = self.dfu.get_objects({'object_refs': [params.get('input_obj_ref')]})
        self._make_fasta(obj_ref=obj['data'][0]['data']['amplicon_set_ref'])
        amp_data = obj['data'][0]['data']

        row_ids = amp_data['data']['row_ids']
        col_ids = amp_data['data']['col_ids']
        values = amp_data['data']['values']
        # Add 'taxonomy' column
        col_ids.append('taxonomy')
        # Make pandas DataFrame
        df = pd.DataFrame(index=row_ids, columns=col_ids)
        for i in range(len(row_ids)):
            df.iloc[i, :-1] = values[i]

        # Get object
        test_row_attributes_permanent_id = obj['data'][0]['data']['row_attributemapping_ref']
        obj = self.dfu.get_objects({'object_refs': [test_row_attributes_permanent_id]})
        tax_dict = obj['data'][0]['data']['instances']

        # Add taxonomy data and transpose matrix
        for row_indx in df.index:
            df.loc[row_indx]['taxonomy'] = tax_dict[row_indx][0]
        return df

    def _get_mdf(self, params):
        """
        Get metadata object and make pd.DataFrame from with samples as index and specified subsetting column
        """

        logging.info('Getting MetadataObject')

        subsetting_field = params.get('subset_field')
        subsetting_field = subsetting_field['meta_group'][0]
        params['subset_field'] = subsetting_field
        # Get object
        obj = self.dfu.get_objects({'object_refs': [params.get('attribute_mapping_obj_ref')]})
        meta_dict = obj['data'][0]['data']['instances']
        attr_l = obj['data'][0]['data']['attributes']

        # Find index of specified category name
        indx = 0
        for i in range(len(attr_l)):
            if attr_l[i]['attribute'] == subsetting_field:
                indx = i
                break
        # Set metadata_samples
        metadata_samples = meta_dict.keys()
        # Make pandas DataFrame
        mdf = pd.DataFrame(index=metadata_samples, columns=[subsetting_field])
        i = 0
        for key, val in meta_dict.items():
            mdf.iloc[i] = val[indx]
            i += 1
        return mdf

    def insert_newlines(self, string, every):
        return '\n'.join(string[i:i + every] for i in range(0, len(string), every))

    def _make_fasta(self, obj_ref):

        logging.info('Making fasta file from AmpliconSet obj: {}'.format(obj_ref))

        set_obj = self.dfu.get_objects({'object_refs': [obj_ref]})
        OTUs = set_obj['data'][0]['data']['amplicons'].keys()
        with open(os.path.join(self.files_folder, "amp_set.fa"), 'w') as fa_file:

            logging.info('Writing to amp_set.fa file')

            for key in OTUs:
                con_str = '>' + key + '\n'
                con_str += self.insert_newlines(set_obj['data'][0]['data']['amplicons'][key]['consensus_sequence'], 60)
                con_str += '\n'
                fa_file.write(con_str)

    def _make_group_dict(self, mdf, subset_field):
        """
        Make dictionary with a subsetting column value as key and samples of that subsetting column value and values
        """

        logging.info('Making grouping dictionary')

        group_dict = {}
        for sample, group in zip(mdf.index, mdf[subset_field]):
            try:
                group_dict[group].append(sample)
            except KeyError:
                group_dict.update({group: [sample]})

        for group, sample_list in group_dict.items():
            group_dict[group].append('taxonomy')

        return group_dict

    def _create_subset_matrices(self, df, mdf, subset_field):
        """
        create dictionary of subset pd.DataFrames
        """

        logging.info('Creating matrices...')

        group_dict = self._make_group_dict(mdf=mdf, subset_field=subset_field)

        # Create dict of sub matrices
        dict_of_sub_matrices = {}
        for key, val in group_dict.items():
            data = df[val]
            dict_of_sub_matrices.update({key: data})

        # Drop rows that have all zero counts
        for key, matrix in dict_of_sub_matrices.items():
            to_drop = []
            for indx in matrix.index:
                if all(val == 0 for val in matrix.loc[indx][0:-1]):
                    to_drop.append(indx)
            dict_of_sub_matrices[key] = matrix.drop(to_drop)

        return dict_of_sub_matrices

    def _save_matrices(self, matrices):
        """
        takes a dictionary of pd.matrices and saves the matrices as tab sep csv's, with the name being keys
        """

        logging.info('Saving matrices: {}'.format(matrices.keys()))

        for group, matrix in matrices.items():
            name = group+'.csv'
            matrix.to_csv(os.path.join(self.files_folder, name), sep='\t')

    def _create_html_report(self):
        """
        Create html report of files in zip by walking through output folder
        """

        logging.info('Creating html report..')

        html_str = '<html>'
        html_str += '<h3>Files In Output Zip File:</h3>\n'
        for root, folders, files in os.walk(self.output_dir):
            # Find the image files by their extensions.
            for f in files:
                if re.match('^[a-zA-Z]+.*.(fa|csv)$', f):  # jpeg|jpg|bmp|png|tiff|pdf|ps|
                    html_str += '<p>' + f + '</p>\n'
        html_str += '</html>'

        with open(os.path.join(self.files_folder, "index.html"), 'w') as index_file:
            index_file.write(html_str)

        # have needed files saved to folder before shock
        shock = self.dfu.file_to_shock({'file_path': self.files_folder,
                                        'make_handle': 0,
                                        'pack': 'zip'})
        # list that goes to 'html_links'
        self.html_paths.append({'shock_id': shock['shock_id'],
                                'name': 'index.html',
                                'label': 'Report',
                                'description': "files in zip"})
        # list that goes to 'file_pahts'
        self.file_paths.append(os.path.join(self.files_folder, 'files.zip'))

    def _call_and_create_objects(self, params):

        logging.info('_call_and_create_objects method')

        list_of_matrix_files = []
        groups = []
        for root, folders, files in os.walk(self.files_folder):

            logging.info('Finding files..')

            # Find the image files by their extensions.
            for f in files:
                if re.match('^[a-zA-Z]+.*.(fa)$', f):
                    fa_file = os.path.join(root, f)
                if re.match('^[a-zA-Z]+.*.(csv)$', f):
                    groups.append(f[0:-4])
                    list_of_matrix_files.append(os.path.join(root, f))

        for csv_file_path, group_name in zip(list_of_matrix_files, groups):

            logging.info('Sending data to importer:\n'
                         'csv_file_path: {}\n'
                         'group_name: {}\n'
                         'fa_file: {}'.format(csv_file_path, group_name, fa_file))

            params['obj_type'] = 'AmpliconMatrix'
            params['matrix_name'] = group_name
            params['tsv_fasta'] = {
                'tsv_file_tsv_fasta': csv_file_path,
                'fasta_file_tsv_fasta': fa_file,
                'metadata_keys_tsv_fasta': 'taxonomy_id, taxonomy, taxonomy_source, consensus_sequence'
            }
            params['scale'] = 'raw'
            params['description'] = 'dsc'
            params['amplicon_set_name'] = group_name+'-set'
            params['sample_set_ref'] = params.get('attribute_mapping_obj_ref')
            params['input_local_file'] = True

            logging.info('Sending params: {}'.format(json.dumps(params, indent=1)))

            obj_run = self.GenAPI.import_matrix_from_biom(params=params)
            logging.info('Object run: {}'.format(obj_run))

    def _create_amp(self):

        amp_structure = {'data': [{'data': {'amplicon_set_ref': '',
                                            'col_attributemapping_ref': '',
                                            'col_mapping': {},
                                            'data': {},
                                            'row_attributemapping_ref': '',
                                            'row_mapping': {},
                                            'scale': 'raw'},
                                  'info': [

                                  ],
                                  'path': [''],
                                  'provenance': [],
                                  'creator': '',
                                  'orig_wsid': 0000,
                                  'created': '',
                                  'epoch': 0000,
                                  'refs': [],
                                  'copy_source_inaccessible': 0,
                                  'extracted_ids': {}}]}

    def run(self, params):

        logging.info('--->\nrunning Amp_Subset_Util with input \n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        df = self._get_df(params)
        mdf = self._get_mdf(params)
        matrices = self._create_subset_matrices(df=df, mdf=mdf, subset_field=params.get('subset_field'))
        self._save_matrices(matrices)
        self._create_html_report()
        self._call_and_create_objects(params)

        return {
            'file_paths': self.file_paths,
            'html_paths': self.html_paths
        }


"""
    arguments:
    obj_type: one of ExpressionMatrix, FitnessMatrix, DifferentialExpressionMatrix
    matrix_name: matrix object name
    workspace_name: workspace name matrix object to be saved to
    input_shock_id: file shock id
    or
    input_file_path: absolute file path
    or
    input_staging_file_path: staging area file path

    optional arguments:
    col_attributemapping_ref: column AttributeMapping reference
    row_attributemapping_ref: row AttributeMapping reference
    genome_ref: genome reference
    matrix_obj_ref: Matrix reference
"""

'''   
typedef structure {
      string obj_type;
      string input_shock_id;
      string input_file_path;
      string input_staging_file_path;
      string matrix_name;
      string amplicon_set_name;
      string scale;
      string description;
      workspace_name workspace_name;

      obj_ref genome_ref;
      obj_ref col_attributemapping_ref;
      obj_ref row_attributemapping_ref;
      obj_ref diff_expr_matrix_ref;
      obj_ref biochemistry_ref;
      obj_ref reads_set_ref;
      obj_ref sample_set_ref;
  } ImportMatrixParams;
'''