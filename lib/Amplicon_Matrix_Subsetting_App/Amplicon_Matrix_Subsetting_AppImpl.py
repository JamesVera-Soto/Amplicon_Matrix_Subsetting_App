# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from Amplicon_Matrix_Subsetting_App.Amp_Subset_Util import Subsetting_Matrices

from installed_clients.KBaseReportClient import KBaseReport
#END_HEADER


class Amplicon_Matrix_Subsetting_App:
    '''
    Module Name:
    Amplicon_Matrix_Subsetting_App

    Module Description:
    A KBase module: Amplicon_Matrix_Subsetting_App
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.config['SDK_CALLBACK_URL'] = self.callback_url
        self.config['KB_AUTH_TOKEN'] = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        self.subsetting_matrices = Subsetting_Matrices(self.config)

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_Amplicon_Matrix_Subsetting_App(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_Amplicon_Matrix_Subsetting_App
        AMS_output = self.subsetting_matrices.run(params)
        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created': [],
                                                'text_message': params['parameter_1']},
                                                'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_Amplicon_Matrix_Subsetting_App

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_Amplicon_Matrix_Subsetting_App return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
