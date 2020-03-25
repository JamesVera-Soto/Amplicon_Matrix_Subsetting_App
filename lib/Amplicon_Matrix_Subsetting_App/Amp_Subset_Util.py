import pandas as pd
import os
import uuid
import logging
import json
import zipfile
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.DataFileUtilClient import DataFileUtil

# Subsetting Class
class Subsetting_Matrices:

    SUBSET_OUT_DIR = 'subsetting_output'

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.scratch = config['scratch']

        self.dfu = DataFileUtil(self.callback_url)

        # set up directory in scratch
        self.output_dir = os.path.join(self.scratch, self.SUBSET_OUT_DIR)
        try:
            os.mkdir(self.output_dir)
        except FileExistsError:
            pass
        # set up directory for files folder
        self.files_folder = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(self.files_folder)

        self.file_paths = []

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
        with open(os.path.join(self.output_dir, "amp_set.fa"), 'w') as fa_file:

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

        return group_dict

    def _create_subset_matrices(self, df, mdf, subset_field):
        """
        create dictionary of subset pd.DataFrames
        """

        logging.info('Creating matrices...')

        group_dict = self._make_group_dict(mdf=mdf, subset_field=subset_field)

        dict_of_sub_matrices = {}
        for key, val in group_dict.items():
            data = df[val]
            dict_of_sub_matrices.update({key: data})

        return dict_of_sub_matrices

    def _save_matrices(self, matrices):
        """
        takes a dictionary of matrices and saves the matrices as tab sep csv's, with the name being keys
        """

        logging.info('Saving matrices: {}'.format(matrices.keys()))

        for group, matrix in matrices.items():
            name = group+'.csv'
            matrix.to_csv(os.path.join(self.output_dir, name), sep='\t')

        self._zip_folder(self.output_dir, os.path.join(self.files_folder, 'Subsetting_output.zip'))

        self.file_paths.append(os.path.join(self.files_folder, 'Subsetting_output.zip'))

    def _zip_folder(self, folder_path, output_path):
        """
        _zip_folder: Zip the contents of an entire folder (with that folder included in the
         archive). Empty subfolders could be included in the archive as well if the 'Included
         all subfolders, including empty ones' portion.
         portion is used.
        """
        with zipfile.ZipFile(output_path, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as ziph:
            for root, folders, files in os.walk(folder_path):
                # Include all subfolders, including empty ones.
                for folder_name in folders:
                    absolute_fpath = os.path.join(root, folder_name)
                    relative_fpath = os.path.join(os.path.basename(root), folder_name)
                    logging.info("Adding folder {} to archive.".format(absolute_fpath))
                    ziph.write(absolute_fpath, relative_fpath)
                for f in files:
                    absolute_path = os.path.join(root, f)
                    relative_path = os.path.join(os.path.basename(root), f)
                    logging.info("Adding file {} to archive.".format(absolute_path))
                    ziph.write(absolute_path, relative_path)

        logging.info("{} created successfully.".format(output_path))

    def run(self, params):

        logging.info('--->\nrunning Amp_Subset_Util with input \n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        df = self._get_df(params)
        mdf = self._get_mdf(params)
        matrices = self._create_subset_matrices(df=df, mdf=mdf, subset_field=params.get('subset_field'))
        self._save_matrices(matrices)
        return {
            'file_paths': self.file_paths
        }
