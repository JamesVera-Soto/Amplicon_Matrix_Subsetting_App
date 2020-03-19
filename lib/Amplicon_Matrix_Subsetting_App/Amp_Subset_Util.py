import pandas as pd
import os
import uuid
import logging
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.DataFIleUtilClient import DataFileUtil

# Subsetting Class
class Subsetting_Matrices:

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.scratch = config['scratch']

        self.dfu = DataFileUtil(self.callback_url)

    def _get_df(self, params):
        logging.info('Getting MatrixObject')
        # Amplicon data
        obj = self.dfu.get_objects({'object_refs': params.get('input_obj_ref')})
        amp_data = obj['data'][0]['data']

        row_ids = amp_data['data']['row_ids']
        col_ids = amp_data['data']['col_ids']
        values = amp_data['data']['values']

        # Make pandas DataFrame
        df = pd.DataFrame(values, index=row_ids, columns=col_ids)
        df = df.T
        return df

    def _get_mdf(self, params):
        logging.info('Getting MetadataObject')
        subsetting_field = params.get('subsetting_field')
        subsetting_field = subsetting_field['meta_group'][0]
        # Get object
        obj = self.dfu.get_objects({'object_refs': params.get('attribute_mapping_obj_ref')})
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

    def create_subset_matrices(self, df, mdf):
        pass

    def run(self, params):
        df = self._get_df(params)
        mdf = self._get_mdf(params)
        self.create_subset_matrices(df=df, mdf=mdf)
