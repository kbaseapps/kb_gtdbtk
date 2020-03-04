# -*- coding: utf-8 -*-
import os
import time
import unittest
from pathlib import Path
from configparser import ConfigParser
from shutil import copyfile

from kb_gtdbtk.kb_gtdbtkImpl import kb_gtdbtk
from kb_gtdbtk.kb_gtdbtkServer import MethodContext
from kb_gtdbtk.authclient import KBaseAuth as _KBaseAuth
from installed_clients.AssemblyUtilClient import AssemblyUtil

from installed_clients.WorkspaceClient import Workspace

# TODO add manual integration tests that *are not* run as part of the test suite
# Full run requires 100MB and 30m at least, infeasable to be part of automated tests
# Clearly document test setup.
# Full test should use kb-sdk test if possible, maybe config option?


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
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL, token=cls.token)
        cls.serviceImpl = kb_gtdbtk(cls.cfg)
        cls.scratch = Path(cls.cfg['scratch']).absolute()
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_gktb_tk_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsid = ret[0]

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_your_method(self):
        tempdir = self.scratch / 'tempstuff'
        tempdir.mkdir(parents=True, exist_ok=True)
        assyfile = tempdir / 'tiny_genome.fa'
        copyfile(Path(__file__).parent / 'tiny_genome.fa', assyfile)

        au = AssemblyUtil(self.callback_url, token=self.token)
        assy = au.save_assembly_from_fasta(
            {'file': {'path': str(assyfile),
                      'assembly_name': 'tiny_genome.fa'
                      },
             'workspace_name': self.wsName,  # TODO AU should take an ID
             'assembly_name': 'tiny_genome.fa'  # well this is redundant
             })

        report = self.serviceImpl.run_kb_gtdbtk(self.ctx, {
            'workspace_id': self.wsid,
            'input_object_ref': assy})
        
        del report
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        # print("START TEST 1\n")
        # ret = self.serviceImpl.run_kb_gtdbtk(self.ctx, {'workspace_name': self.wsName,
        #                                                      'inputObjectRef': '27079/17/1'})
