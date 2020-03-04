# -*- coding: utf-8 -*-
import hashlib
import os
import requests
import time
import unittest
import uuid
from pathlib import Path
from configparser import ConfigParser
from shutil import copyfile

from kb_gtdbtk.kb_gtdbtkImpl import kb_gtdbtk
from kb_gtdbtk.kb_gtdbtkServer import MethodContext
from kb_gtdbtk.authclient import KBaseAuth as _KBaseAuth
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.AbstractHandleClient import AbstractHandle

# TODO add manual integration tests that *are not* run as part of the test suite
# Full run requires 100MB and 30m at least, infeasable to be part of automated tests
# Clearly document test setup.
# Full test should use kb-sdk test if possible, maybe config option?


# TODO add tests for other types
# TODO add tests for copy of a copy
# TODO add a few failing case tests (most covered by unit tests)
# each test takes 30m ... need to think of a way to speed this up. Maybe a flag to skip running
# GTDB-tk so the report is empty and just check the downloads?

WORKDIR = '/kb/module/work/tmp/'


class kb_gtdbtkTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_gtdbtk'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_gtdbtk',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.ws_url = cls.cfg['workspace-url']
        cls.shock_url = cls.cfg['shock-url']
        cls.ws = Workspace(cls.ws_url, token=cls.token)
        cls.serviceImpl = kb_gtdbtk(cls.cfg)
        cls.scratch = Path(cls.cfg['scratch']).absolute()
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_gktb_tk_" + str(suffix)
        ret = cls.ws.create_workspace({'workspace': cls.wsName})
        cls.wsid = ret[0]
        cls.hs = AbstractHandle(cls.cfg['handle-service-url'], token=cls.token)
        cls.au = AssemblyUtil(cls.callback_url, token=cls.token)
        cls.dfu = DataFileUtil(cls.callback_url, token=cls.token)
        cls.handles_to_delete = []
        cls.nodes_to_delete = []

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.ws.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'nodes_to_delete'):
            for node in cls.nodes_to_delete:
                cls.delete_shock_node(node)
        if hasattr(cls, 'handles_to_delete'):
            if cls.handles_to_delete:
                cls.hs.delete_handles(cls.hs.hids_to_handles(cls.handles_to_delete))
                print('Deleted handles ' + str(cls.handles_to_delete))

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.token)}
        requests.delete(cls.shock_url + '/node/' + node_id, headers=header, allow_redirects=True)
        print('Deleted shock node ' + node_id)

    def test_your_method(self):
        tempdir = self.scratch / 'tempstuff'
        tempdir.mkdir(parents=True, exist_ok=True)
        assyfile = tempdir / 'tiny_genome.fa'
        copyfile(Path(__file__).parent / 'tiny_genome.fa', assyfile)

        assy = self.au.save_assembly_from_fasta(
            {'file': {'path': str(assyfile)},
             'workspace_name': self.wsName,  # TODO AU should take an ID
             'assembly_name': 'tiny_genome.fa'
             })

        report = self.serviceImpl.run_kb_gtdbtk(self.ctx, {
            'workspace_id': self.wsid,
            'input_object_ref': assy})[0]

        md5s = {
            'index.html': 'e865d72e375bbbc5721f8d999698e1c5'
        }

        self.check_gtdbtk_output(report, 2277, md5s)

    def check_gtdbtk_output(self, report, zipsize, filename_to_md5, tolerance=15):
        print(report)
        minsize = zipsize - tolerance
        maxsize = zipsize + tolerance
        filename = 'output.zip'

        objname = report['report_name']
        obj = self.dfu.get_objects({'object_refs': [report['report_ref']]})['data'][0]

        print(obj)
        d = obj['data']
        shock_url = d['html_links'][0]['URL']
        hid = d['html_links'][0]['handle']
        del d['html_links'][0]['URL']
        del d['html_links'][0]['handle']

        assert objname == obj['info'][1]
        assert d == {'direct_html': None,
                     'direct_html_link_index': 0,
                     'file_links': [],
                     'html_links': [{'description': 'HTML report for GTDBTk Classify',
                                     'label': 'index.html',
                                     'name': 'index.html'}],
                     'html_window_height': None,
                     'objects_created': [],
                     'summary_window_height': None,
                     'text_message': None,
                     'warnings': []}

        shocknode = shock_url.split('/')[-1]
        self.handles_to_delete.append(hid)
        self.nodes_to_delete.append(shocknode)

        shockret = requests.get(self.shock_url + '/node/' + shocknode,
                                headers={'Authorization': 'OAuth ' + self.token}).json()['data']
        print(shockret)

        assert shockret['id'] == shocknode
        shockfile = shockret['file']
        assert shockfile['name'] == filename
        assert shockfile['size'] > minsize  # zip file size & md5 not repeatable
        assert shockfile['size'] < maxsize

        handleret = self.hs.hids_to_handles([hid])[0]
        print(handleret)
        assert handleret['url'] == self.shock_url
        assert handleret['hid'] == hid
        assert handleret['file_name'] == filename
        assert handleret['type'] == 'shock'
        assert handleret['id'] == shocknode
        assert handleret['remote_md5'] == shockfile['checksum']['md5']

        # check data in shock
        zipdir = Path(WORKDIR) / str(uuid.uuid4())
        self.dfu.shock_to_file(
            {'shock_id': shocknode,
             'unpack': 'unpack',
             'file_path': str(zipdir / filename)
             })
        files = os.listdir(zipdir)
        files.remove(filename)
        print(files)

        assert files == list(filename_to_md5.keys())
        for f in files:
            assert hashlib.md5(open(zipdir / f, 'rb').read()).hexdigest() == filename_to_md5[f]
