from pytest import raises

from unittest.mock import create_autospec
from pathlib import Path

from kb_gtdbtk.core.kb_client_set import KBClients
from kb_gtdbtk.core.sequence_downloader import download_sequence
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.WorkspaceClient import Workspace
from core.test_utils import assert_exception_correct


def _client_mocks():
    clis = create_autospec(KBClients, spec_set=True, instance=True)
    dfu = create_autospec(DataFileUtil, spec_set=True, instance=True)
    au = create_autospec(AssemblyUtil, spec_set=True, instance=True)
    mgu = create_autospec(MetagenomeUtils, spec_set=True, instance=True)
    ws = create_autospec(Workspace, spec_set=True, instance=True)

    clis.dfu.return_value = dfu
    clis.au.return_value = au
    clis.mgu.return_value = mgu
    clis.ws.return_value = ws
    return clis


def test_with_sets_genomeset():
    _check_genomeset('KBaseSets.GenomeSet-1.0', {'items': [{'ref': '1/2/3'}, {'ref': '4/5/6'}]})


def test_with_search_genomeset():
    _check_genomeset(
        'KBaseSearch.GenomeSet-2.1',
        # The search typespec doesn't describe the semantics of the inner map keys
        {'elements': {'some_key': {'ref': '1/2/3'}, 'another_key': {'ref': '4/5/6'}}})


def _check_genomeset(type_, data):
    clis = _client_mocks()

    clis.dfu().get_objects.return_value = {
        'data': [
            {'info': [
                3,
                'myobj',
                type_,
                'ignored time',
                7,
                'someuser',
                34567,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': data
             }
        ]
    }

    clis.ws().get_objects2.side_effect = [
        {'data': [
            {'info': [
                2,
                'myobj',
                'KBaseGenomes.Genome-11.0',
                'ignored time',
                1,
                'someuser',
                3,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': {'contigset_ref': '45/21/89'}
             }
            ]
         },
        {'data': [
            {'info': [
                5,
                'myobj',
                'KBaseGenomes.Genome-11.0',
                'ignored time',
                6,
                'someuser',
                4,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': {'assembly_ref': '7/8/9'}
             }
            ]
         },
        None  # cause fail if called too many times
    ]

    clis.au().get_assembly_as_fasta.side_effect = [
        {'path': '/path/to/file1', 'assembly_name': 'assyname1'},
        {'path': '/path/to/file2', 'assembly_name': 'assyname2'},
        None  # cause fail if called too many times
    ]

    ret = download_sequence('34567/3/7', Path('somepath/or/other'), clis)

    assert ret == {'45_21_89': {'path': '/path/to/file1', 'assembly_name': 'assyname1'},
                   '7_8_9': {'path': '/path/to/file2', 'assembly_name': 'assyname2'}
                   }

    assert clis.dfu().get_objects.call_args_list == [(({'object_refs': ['34567/3/7']},), {})]

    assert clis.ws().get_objects2.call_args_list == [
        (({'objects': [{'ref': '34567/3/7;1/2/3'}]},), {}),
        (({'objects': [{'ref': '34567/3/7;4/5/6'}]},), {})]

    assert clis.au().get_assembly_as_fasta.call_args_list == [
        (({'ref': '34567/3/7;1/2/3;45/21/89', 'filename': 'somepath/or/other/45_21_89'},), {}),
        (({'ref': '34567/3/7;4/5/6;7/8/9', 'filename': 'somepath/or/other/7_8_9'},), {})]


def test_with_genome():
    _check_genome('assembly_ref')
    _check_genome('contigset_ref')


def _check_genome(ref_key):
    clis = _client_mocks()

    obj_data = {
        'data': [
            {'info': [
                3,
                'myobj',
                'KBaseGenomes.Genome-11.0',
                'ignored time',
                7,
                'someuser',
                34567,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': {ref_key: '1/2/3'}
             }
        ]
    }

    clis.dfu().get_objects.return_value = obj_data

    # duplicate call, could be optimized out with increase in code complexity
    clis.ws().get_objects2.return_value = obj_data

    clis.au().get_assembly_as_fasta.return_value = {
        'path': '/path/to/file3', 'assembly_name': 'assyname3'}

    ret = download_sequence('34567/3/7', Path('somepath/or/other'), clis)

    assert ret == {'1_2_3': {'path': '/path/to/file3', 'assembly_name': 'assyname3'}
                   }

    assert clis.dfu().get_objects.call_args_list == [(({'object_refs': ['34567/3/7']},), {})]

    assert clis.ws().get_objects2.call_args_list == [(({'objects': [{'ref': '34567/3/7'}]},), {})]

    assert clis.au().get_assembly_as_fasta.call_args_list == [
        (({'ref': '34567/3/7;1/2/3', 'filename': 'somepath/or/other/1_2_3'},), {})]


def test_with_assembly():
    _check_assembly('KBaseGenomes.ContigSet-3.1')
    _check_assembly('KBaseGenomeAnnotations.Assembly-5.0')


def _check_assembly(type_):
    clis = _client_mocks()

    obj_data = {
        'data': [
            {'info': [
                3,
                'myobj',
                type_,
                'ignored time',
                7,
                'someuser',
                34567,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': {}  # no data needed, could optimize by just calling get info vs get object
             }
        ]
    }

    clis.dfu().get_objects.return_value = obj_data

    clis.au().get_assembly_as_fasta.return_value = {
        'path': '/path/to/file4', 'assembly_name': 'assyname4'}

    ret = download_sequence('34567/3/7', Path('somepath/or/other'), clis)

    assert ret == {'34567_3_7': {'path': '/path/to/file4', 'assembly_name': 'assyname4'}
                   }

    assert clis.dfu().get_objects.call_args_list == [(({'object_refs': ['34567/3/7']},), {})]

    assert clis.au().get_assembly_as_fasta.call_args_list == [
        (({'ref': '34567/3/7', 'filename': 'somepath/or/other/34567_3_7'},), {})]


def test_with_assemblyset():
    clis = _client_mocks()

    clis.dfu().get_objects.return_value = {
        'data': [
            {'info': [
                3,
                'myobj',
                'KBaseSets.AssemblySet-1.0',
                'ignored time',
                7,
                'someuser',
                34567,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': {'items': [{'ref': '1/2/3'}, {'ref': '4/5/6'}]}
             }
        ]
    }

    clis.ws().get_objects2.side_effect = [
        {'data': [
            {'info': [
                2,
                'myobj',
                'KBaseGenomeAnnotations.Assembly-5.0',
                'ignored time',
                1,
                'someuser',
                3,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': {'contigset_ref': '45/21/89'}
             }
            ]
         },
        {'data': [
            {'info': [
                5,
                'myobj',
                'KBaseGenomeAnnotations.Assembly-5.0',
                'ignored time',
                6,
                'someuser',
                4,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': {'assembly_ref': '7/8/9'}
             }
            ]
         },
        None  # cause fail if called too many times
    ]

    clis.au().get_assembly_as_fasta.side_effect = [
        {'path': '/path/to/file1', 'assembly_name': 'assyname1'},
        {'path': '/path/to/file2', 'assembly_name': 'assyname2'},
        None  # cause fail if called too many times
    ]

    ret = download_sequence('34567/3/7', Path('somepath/or/other'), clis)

    assert ret == {'1_2_3': {'path': '/path/to/file1', 'assembly_name': 'assyname1'},
                   '4_5_6': {'path': '/path/to/file2', 'assembly_name': 'assyname2'}
                   }

    assert clis.dfu().get_objects.call_args_list == [(({'object_refs': ['34567/3/7']},), {})]

    assert clis.au().get_assembly_as_fasta.call_args_list == [
        (({'ref': '34567/3/7;1/2/3', 'filename': 'somepath/or/other/1_2_3'},), {}),
        (({'ref': '34567/3/7;4/5/6', 'filename': 'somepath/or/other/4_5_6'},), {})]


def test_fail_bad_type():
    clis = _client_mocks()

    clis.dfu().get_objects.return_value = {
        'data': [
            {'info': [
                3,
                'myobj',
                'KBaseGenomes.GenomeComparison-2.1',
                'ignored time',
                7,
                'someuser',
                34567,
                'some_workspace',
                'md5 here',
                416467216,
                {}
                ],
             'data': {'items': [{'ref': '1/2/3'}, {'ref': '4/5/6'}]}
             }
        ]
    }

    with raises(Exception) as got:
        download_sequence('34567/3/7', Path('somepath/or/other'), clis)
    assert_exception_correct(
        got.value, ValueError('KBaseGenomes.GenomeComparison type is not supported'))
