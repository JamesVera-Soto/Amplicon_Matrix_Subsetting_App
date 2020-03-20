/*
A KBase module: Amplicon_Matrix_Subsetting_App
*/

module Amplicon_Matrix_Subsetting_App {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    typedef structure {
        string input_obj_ref;
        string attribute_mapping_obj_ref;
        mapping<string, string> subset_field;
    } Amp_Subset_Params;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_Amplicon_Matrix_Subsetting_App(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
