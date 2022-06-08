import uuid
import tempfile

from pathlib import Path
from unittest.mock import create_autospec

from kb_gtdbtk.core.kb_report_generation import generate_report
from kb_gtdbtk.core.kb_client_set import KBClients
from installed_clients.KBaseReportClient import KBaseReport


def _client_mocks():
    clis = create_autospec(KBClients, spec_set=True, instance=True)
    report = create_autospec(KBaseReport, spec_set=True, instance=True)

    clis.report.return_value = report
    return clis


def test_create_report():
    clis = _client_mocks()

    with tempfile.TemporaryDirectory(prefix='test_create_report') as test_dir_str:
        test_dir = Path(test_dir_str)

        clis.report().create_extended_report.return_value = {
            'ref': '78/5/1',
            'name': 'GTDBTk_report_bd5c9ba0-db2c-4c03-b2de-92d9c94ce51e'
            }

        ret = generate_report(
            clis, test_dir, 78, None,
            uuid_gen=lambda: uuid.UUID('bd5c9ba0-db2c-4c03-b2de-92d9c94ce51e'))

        shock_id = ret['archive_shock_id']
        
        assert ret == {
            'report_ref': '78/5/1',
            'report_name': 'GTDBTk_report_bd5c9ba0-db2c-4c03-b2de-92d9c94ce51e',
            'archive_shock_id': shock_id
            }

        clis.report().create_extended_report.assert_called_once_with(
                {'direct_html_link_index': 0,
                 'html_links': [{'path': str(test_dir),
                                 'name': 'index.html',
                                 'label': 'index.html',
                                 'description': 'HTML report for GTDBTk Classify'
                                 }],
                 'file_links': [{'shock_id': shock_id,
                                 'name': 'GTDB-Tk_classify_wf.zip',
                                 'description': 'GTDB-Tk Classify WF output'
                                 }],
                 'report_object_name': 'GTDBTk_report_bd5c9ba0-db2c-4c03-b2de-92d9c94ce51e',
                 'workspace_id': 78,
                 })

        assert (test_dir / 'index.html').is_file()
