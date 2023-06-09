from pytest import raises
from kb_gtdbtk.core.api_translation import get_gtdbtk_params
from core.test_utils import assert_exception_correct


def test_get_gtdbtk_params():
    p = get_gtdbtk_params({
        'workspace_id': 56,
        'input_object_ref': '5/6/7',
        'inputObjectRef': 'should be ignored',  # old key name
        'copy_proximals': 0,
        'some random key': 'foo'  # should this be an error?
    })
    assert p == ('5/6/7', 56, 0, 10, 0, 0, 0)

    p = get_gtdbtk_params({'workspace_id': 92, 'input_object_ref': '104/67/3', 'min_perc_aa': 78.9, 'copy_proximals': 0})
    assert p == ('104/67/3', 92, 0, 78.9, 0, 0, 0)


def test_get_gtdbtk_params_backwards_compatibility():
    p = get_gtdbtk_params({'workspace_id': 56, 'inputObjectRef': '8/9/10', 'copy_proximals': 0})
    assert p == ('8/9/10', 56, 0, 10, 0, 0, 0)


def test_get_gtdbtk_params_fail_bad_args():
    _get_gtdbtk_params_fail({'workspace_id': 7}, ValueError(
        'input_object_ref is required and must be a string'))
    _get_gtdbtk_params_fail({'input_object_ref': 78, 'workspace_id': 7}, ValueError(
        'input_object_ref is required and must be a string'))
    _get_gtdbtk_params_fail({'inputObjectRef': 78, 'workspace_id': 7}, ValueError(
        'input_object_ref is required and must be a string'))
    _get_gtdbtk_params_fail(
        {'inputObjectRef': '1/1/1', 'min_perc_aa': 'foo', 'workspace_id': 7},
        ValueError('min_perc_aa must be a float'))
    _get_gtdbtk_params_fail({'inputObjectRef': '1/1/1'}, ValueError(
        'workspace_id is required and must be an integer > 0'))
    _get_gtdbtk_params_fail({'inputObjectRef': '1/1/1', 'workspace_id': 'foo'}, ValueError(
        'workspace_id is required and must be an integer > 0'))
    _get_gtdbtk_params_fail({'inputObjectRef': '1/1/1', 'workspace_id': 0}, ValueError(
        'workspace_id is required and must be an integer > 0'))


def _get_gtdbtk_params_fail(params, expected):
    with raises(Exception) as got:
        get_gtdbtk_params(params)
    assert_exception_correct(got.value, expected)
