/*
A KBase module: kb_gtdbtk
*/

module kb_gtdbtk {

    /* Parameters for the GTDB-tk Classify (classify_wf) run.

        Required:
        input_object_ref: A reference to the workspace object to process.
        workspace_id: The integer workspace ID where the results will be saved.

        Optional:
        min_perc_aa: the minimum sequence alignment as a percent, default 10.
        
    */
    typedef structure {
        string input_object_ref;
        int workspace_id;
        float min_perc_aa;
    } GTDBtk_Classify_Params;

    
    /* The results of the GTDB-tk run.

        report_name: The name of the report object in the workspace.
        report_ref: The UPA of the report object, e.g. wsid/objid/ver.
    */
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;


    /*
        Run GTDB-tk Classify
    */
    funcdef run_kb_gtdbtk_classify_wf(GTDBtk_Classify_Params params) returns (ReportResults output)
        authentication required;

};
