'''
Functions to translate values between the SDK API and core app logic.
'''

from typing import Dict, NamedTuple as _NamedTuple

class GTDBTKParams(_NamedTuple):
    '''
    Parameters for running GTDB-tk in the KBase environment.
    '''

    ref: str
    ''' The workspace object reference in the x/y/z form. '''

    workspace_id: int
    ''' The ID of the KBase workspace where the report should be saved. '''

    min_perc_aa: int
    ''' The mimimum sequence alignment in percent. '''
    

def get_gtdbtk_params(input_params: Dict[str, object]) -> GTDBTKParams:
    '''
    Process input parameters supplied to the GTDB-tk run method and parse them into
    a named tuple. The expected fields are:

    input_object_ref: a workspace reference to the input object that will be processed in the
        x/y/z form.
    workspace_id: the integer ID of the workspace where the results will be saved.
    min_perc_aa: the minimum sequence alignment as a percent.

    :param input_params: the input parameters passed to the method.
    :returns: the parsed parameters.
    :raises ValueError: if any of the parameters are invalid.
    '''
    # TODO put actual params in spec
    ref = input_params.get('input_object_ref')  # TODO Update UI
    if not ref:
        # for backwards compatibility
        ref = input_params.get('inputObjectRef')
    if type(ref) != str:
        raise ValueError('input_object_ref is required and must be a string')
    # could check ref format, but the ws will do that for us. YAGNI.
    
    min_perc_aa = input_params.get('min_perc_aa', 10)
    if type(min_perc_aa) != int:
        raise ValueError('min_perc_aa must be an integer')
    # TODO check 0 <= min_perc_aa <= 1
    
    wsid = input_params.get('workspace_id')
    if type(wsid) != int or wsid < 1:
        raise ValueError('workspace_id is required and must be an integer > 0')

    return GTDBTKParams(ref, wsid, min_perc_aa)
