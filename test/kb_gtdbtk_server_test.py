# -*- coding: utf-8 -*-
import hashlib
import os
import requests
import time
import unittest
import uuid
import json
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

# TODO add tests for other types (and make some of them actually find genes)
# TODO add tests for copy of a copy
# TODO add a few failing case tests (most covered by unit tests)

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
        cls.prepare_data()

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

    @classmethod
    def ref_from_info(cls, obj_info):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        return "/".join([str(obj_info[WSID_I]), str(obj_info[OBJID_I]), str(obj_info[VERSION_I])])
        
    @classmethod
    def prepare_data(cls):
        tempdir = cls.scratch / 'tempstuff'
        tempdir.mkdir(parents=True, exist_ok=True)

        # single assembly
        this_filename = 'Rhodo_contigs.fa'
        single_assyfile = tempdir / this_filename
        copyfile(Path(__file__).parent / 'data' / this_filename, single_assyfile)
        cls.single_assy = cls.au.save_assembly_from_fasta(
            {'file': {'path': str(single_assyfile)},
             'workspace_name': cls.wsName,  # TODO AU should take an ID
             'assembly_name': this_filename
             })

        # MG assembly
        this_filename = '37AB_metaSPAdes_binnedcontigs.contigs.gz'
        mg_assyfile = tempdir / this_filename
        copyfile(Path(__file__).parent / 'data' / this_filename, mg_assyfile)
        mg_assy = cls.au.save_assembly_from_fasta(
            {'file': {'path': str(mg_assyfile)},
             'workspace_name': cls.wsName,
             'assembly_name': this_filename
             })

        # binned contigs
        this_filename = '37AB_metaSPAdes.binnedcontigs'
        mg_binfile = tempdir / this_filename
        copyfile(Path(__file__).parent / 'data' / this_filename, mg_binfile)
        with open (mg_binfile, 'r') as bin_fh:
            bin_obj = json.load(bin_fh)
        bin_obj['assembly_ref'] = mg_assy
        bin_obj_info = cls.ws.save_objects({
            'workspace': cls.wsName,
            'objects': [
                {
                    'type': 'KBaseMetagenomes.BinnedContigs',
                    'data': bin_obj,
                    'name': this_filename
                }
            ]})[0]
        cls.binned_contigs = cls.ref_from_info(bin_obj_info)
        
        
    # test binnedcontigs input
    # NOT ABLE TO RUN ON DEV1.  Too much memory required
    @unittest.skip("skipped test_classify_wf_binnedcontigs()")  # uncomment to skip
    def test_classify_wf_binnedcontigs(self):
        report = self.serviceImpl.run_kb_gtdbtk_classify_wf(self.ctx, {
            'workspace_id': self.wsid,
            'input_object_ref': self.binned_contigs})[0]
        # TODO: after shrinking data to fit on dev1, test report content
        

    # test assembly input
    # HIDE @unittest.skip("skipped test_classify_wf_assembly()")  # uncomment to skip
    def test_classify_wf_assembly(self):
        report = self.serviceImpl.run_kb_gtdbtk_classify_wf(self.ctx, {
            'workspace_id': self.wsid,
            'input_object_ref': self.single_assy})[0]

        # can't easily maintain md5s through repeated updates.  don't require
        md5s = {}
        """
        # this is for module v0.1.5, GTDB-Tk v1.3.0
        md5s = {
            'index.html': 'e865d72e375bbbc5721f8d999698e1c5',
            'gtdbtk.bac120.markers_summary.tsv': 'a09e7128e6af6d0ff808436a12692777',
            'gtdbtk.bac120.markers_summary.tsv.json': '994dd707e876157f2ad0d9150cb2dc0a',
            'gtdbtk.ar122.markers_summary.tsv': '3d15217ca2e43d27a01327cd7d23a586',
            'gtdbtk.ar122.markers_summary.tsv.json': '124843868858867aba1f43b51707864b',
        }
        """

        self.check_gtdbtk_output(report, 4624, md5s)


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
        # can't maintain filesize expectation through wrapped tool/db updates
        """
        assert shockfile['size'] > minsize  # zip file size & md5 not repeatable
        assert shockfile['size'] < maxsize
        """
        
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

        # can't easily maintain md5s through repeated updates.  don't require
        """
        assert set(files) == set(filename_to_md5.keys())
        for f in files:
            assert hashlib.md5(open(zipdir / f, 'rb').read()).hexdigest() == filename_to_md5[f]
        """
