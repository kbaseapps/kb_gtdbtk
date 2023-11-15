'''
Container for a set of KBase clients for KB-SDK modules. Given a url for the SDK callback server
and a token, initializes the client set.
'''

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.WorkspaceClient import Workspace
from installed_clients.SetAPIClient import SetAPI
from installed_clients.AbstractHandleClient import AbstractHandle

# DEV NOTES: This is not tested in travis and must be tested manually.


class KBClients:
    '''
    A set of clients for KB-SDK modules.
    '''

    def __init__(self, callback_url: str, workspace_url: str, handle_srv_url: str, user_token: str):
        '''
        Create the client set.

        :param callback_url: The url of the callback server.
        :param workspace_url: The url of the KBase workspace server.
        :param user_token: The user's token.
        '''
        # TODO check inputs aren't None or empty string
        self._dfu = DataFileUtil(callback_url, token=user_token)
        self._au = AssemblyUtil(callback_url, token=user_token)
        self._mgu = MetagenomeUtils(callback_url, token=user_token)
        self._report = KBaseReport(callback_url, token=user_token)
        self._ws = Workspace(workspace_url, token=user_token)
        self._setAPI = SetAPI(callback_url, token=user_token)
        self._hs = AbstractHandle(handle_srv_url, token=user_token)

    # Using methods rather than instance variables since create_autospec doesn't play nicely
    # with instance variables.

    def dfu(self):
        '''
        Get the DataFileUtil client.
        :returns: the client.
        '''
        return self._dfu

    def au(self):
        '''
        Get the AssemblyUtil client.
        :returns: the client.
        '''
        return self._au

    def mgu(self):
        '''
        Get the MetagenomUtils client.
        :returns: the client.
        '''
        return self._mgu

    def ws(self):
        '''
        Get the Workspace client.
        :returns: the client.
        '''
        return self._ws

    def report(self):
        '''
        Get the KBaseReport client.
        :returns: the client.
        '''
        return self._report

    def setAPI(self):
        '''
        Get the SetAPI client.
        :returns: the client.
        '''
        return self._setAPI

    def hs(self):
        '''
        Get the HandleService client.
        :returns: the client.
        '''
        return self._hs
