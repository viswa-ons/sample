import os
import sys
import gzip
from datetime import datetime
import time
import logging
import re

from tqdm import tqdm

import numpy as np
import xml.etree.cElementTree as etree

import defaults
import models
from constants import PYUNIPROT_DATA_DIR, PYUNIPROT_DIR

from configparser import RawConfigParser

from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.engine import reflection
from sqlalchemy.sql import sqltypes

if sys.version_info[0] == 3:
    from urllib.request import urlretrieve
    from requests.compat import urlparse, urlsplit
    from subprocess import getoutput
else:
    from urllib import urlretrieve
    from urlparse import urlparse, urlsplit
    from commands import getoutput

log = logging.getLogger(__name__)

alchemy_pandas_dytpe_mapper = {
    sqltypes.Text: np.unicode,
    sqltypes.String: np.unicode,
    sqltypes.Integer: np.float,
    sqltypes.REAL: np.double
}


def get_connection_string(connection=None):
    """return SQLAlchemy connection string if it is set

    :param connection: get the SQLAlchemy connection string #TODO
    :rtype: str
    """
    if not connection:
        config = configparser.ConfigParser()
        cfp = defaults.config_file_path
        if os.path.exists(cfp):
            log.info('fetch database configuration from %s', cfp)
            config.read(cfp)
            connection = config['database']['sqlalchemy_connection_string']
            log.info('load connection string from %s: %s', cfp, connection)
        else:
            with open(cfp, 'w') as config_file:
                connection = defaults.sqlalchemy_connection_string_default
                config['database'] = {'sqlalchemy_connection_string': connection}
                config.write(config_file)
                log.info('create configuration file %s', cfp)

    return connection

class DbManager():
    pmids = set()
    keywords = {}
    subcellular_locations = {}
    tissues = {}

    def __init__(self, connection=None):
        """The DbManager implements all function to upload CTD files into the database. Prefered SQL Alchemy 
        database is MySQL with pymysql.
        
        :param connection: custom database connection SQL Alchemy string
        :type connection: str
        """
        #super(DbManager, self).__init__(connection=connection)
        print("Inside __init__ function")

    def db_import_xml(self, url=None, force_download=False, taxids=None, silent=False):
        """Updates the CTD database
        
        1. downloads gzipped XML
        2. drops all tables in database
        3. creates all tables in database
        4. import XML
        5. close session

        :param taxids: list of NCBI taxonomy identifier
        :type taxids: list
        :param url: iterable of URL strings
        :type url: str
        :param force_download: force method to download
        :type: bool
        """

        log.info('Update CTD database from {}'.format(url))

        #self._drop_tables()
        xml_gzipped_file_path, version_file_path = self.download(url, force_download)
        #self._create_tables()
        #self.import_version(version_file_path)
        self.import_xml(xml_gzipped_file_path, taxids, silent)
        #self.session.close()

    def import_version(self, version_file_path):
        pattern = "UniProtKB/(?P<knowledgebase>Swiss-Prot|TrEMBL) Release" \
                  " (?P<release_name>\\d{4}_\\d{2}) of (?P<release_date>\\d{2}-\\w{3}-\\d{4})"
        with open(version_file_path) as fd:
            content = fd.read()

        for knowledgebase, release_name, release_date_str in re.findall(pattern, content):
            release_date = datetime.strptime(release_date_str, '%d-%b-%Y')

            version = models.Version(
                knowledgebase=knowledgebase,
                release_name=release_name,
                release_date=release_date
            )

            self.session.add(version)

        self.session.commit()

    def import_xml(self, xml_gzipped_file_path, taxids, silent=False):
        """Imports XML

        :param str xml_gzipped_file_path: path to XML file
        :param int,(int,),[int,] taxids: NCBI taxonomy identifier
        :param bool silent: no output if True
        """
        #version = self.session.query(models.Version).filter(models.Version.knowledgebase == 'Swiss-Prot').first()
        #version.import_start_date = datetime.now()

        entry_xml = '<entries>'
        number_of_entries = 0
        interval = 10
        start = False

        if sys.platform in ('linux', 'linux2', 'darwin'):
            log.info('Load gzipped XML from {}'.format(xml_gzipped_file_path))

            zcat_command = 'gzcat' if sys.platform == 'darwin' else 'zcat'

            number_of_lines = int(getoutput("{} {} | wc -l".format(zcat_command, xml_gzipped_file_path)))

            tqdm_desc = 'Import {} lines'.format(number_of_lines)

        else:
            print('bin was anderes')
            number_of_lines = None
            tqdm_desc = None

        xmlns_re = re.compile(' xmlns="[^"]+"')

        with gzip.open(xml_gzipped_file_path) as fd:

            for line in tqdm(fd, desc=tqdm_desc, total=number_of_lines, mininterval=1, disable=silent):

                end_of_file = line.startswith(b"</uniprot>")

                if line.startswith(b"<entry "):
                    start = True

                elif end_of_file:
                    start = False

                #if start:
                #    entry_xml += line.decode("utf-8")

                if start:
                    if line.startswith(b"</entry>") or end_of_file:
                        number_of_entries += 1
                        #start = False

                        if number_of_entries == interval or end_of_file:

                            entry_xml += "</entry></entries>"
                            #print(entry_xml)
                            self.insert_entries(entry_xml, taxids)
                            exit()

                            if end_of_file:
                                break

                            else:
                                entry_xml = "<entries>" + xmlns_re.sub('',line[8:].decode("utf-8"))
                                print(line[8:].decode("utf-8"))
                                number_of_entries = 0
                        else:
                            entry_xml += xmlns_re.sub('',line.decode("utf-8"))
                    else:
                        entry_xml += xmlns_re.sub('',line.decode("utf-8"))

        #version.import_completed_date = datetime.now()
        #self.session.commit()

    def insert_entries(self, entries_xml, taxids):
        """
        insert UniProt entries from XML

        :param str entries_xml: XML string
        :param int,tuple,list taxids: NCBI taxonomy IDs
        :return:
        """

        entries = etree.fromstring(entries_xml)
        del entries_xml

        for entry in entries:
            #print(etree.tostring(entry, 'utf-8', method="xml"))
            self.insert_entry(entry, taxids)
            entry.clear()
            del entry

        entries.clear()
        del entries

        #self.session.commit()

    # profile
    def insert_entry(self, entry, taxids):
        """
        insert UniProt entry"

        :param entry: XML node entry
        :param taxids: int,tuple,list taxids: NCBI taxonomy IDs
        :return:
        """
        entry_dict = entry.attrib
        entry_dict['created'] = datetime.strptime(entry_dict['created'], '%Y-%m-%d')
        entry_dict['modified'] = datetime.strptime(entry_dict['modified'], '%Y-%m-%d')

        #print(entry_dict['created'])

        taxid = self.get_taxid(entry)

        if taxids is None or taxid in taxids:
            #entry_dict = self.update_entry_dict(entry, entry_dict, taxid)
            entry_obj = models.Entry(**entry_dict)
            del entry_dict

            #self.session.add(entry_obj)

    @classmethod
    def get_taxid(cls, entry):
        """
        get NCBI taxonomy identifier from XML node entry

        :param entry:X ML node entry
        :return: int
        """
        query = "./organism/dbReference[@type='NCBI Taxonomy']"
        return int(entry.find(query).get('id'))

    @classmethod
    def download(cls, url=None, force_download=False):
        """Downloads uniprot_sprot.xml.gz and reldate.txt (release date information) from URL or file path

        .. note::

            only URL/path of xml.gz is needed and valid value for parameter url. URL/path for reldate.txt have to be the
            same folder
    
        :param str url: UniProt gzipped URL or file path
        :param force_download: force method to download
        :type force_download: bool
        """
        if url:
            version_url = os.path.join(os.path.dirname(url), defaults.VERSION_FILE_NAME)
        else:
            url = os.path.join(defaults.XML_DIR_NAME, defaults.SWISSPROT_FILE_NAME)
            version_url = os.path.join(defaults.XML_DIR_NAME, defaults.VERSION_FILE_NAME)
        #print(url)
        #print(version_url)
        xml_file_path = cls.get_path_to_file_from_url(url)
        version_file_path = cls.get_path_to_file_from_url(version_url)

        if force_download or not os.path.exists(xml_file_path):

            log.info('download {} and {}'.format(xml_file_path, version_file_path))

            scheme = urlsplit(url).scheme

            if scheme in ('ftp', 'http'):
                urlretrieve(version_url, version_file_path)
                #urlretrieve(url, xml_file_path)

            elif not scheme and os.path.isfile(url):
                shutil.copyfile(url, xml_file_path)
                shutil.copyfile(version_url, version_file_path)

        print(version_file_path)
        print(xml_file_path)
        return xml_file_path, version_file_path

    @classmethod
    def get_path_to_file_from_url(cls, url):
        """standard file path
        
        :param str url: download URL
        """
        file_name = urlparse(url).path.split('/')[-1]
        return os.path.join(PYUNIPROT_DATA_DIR, file_name)

def update(connection=None, urls=None, force_download=False, taxids=None, silent=False):
    """Updates CTD database

    :param urls: list of urls to download
    :type urls: iterable
    :param connection: custom database connection string
    :type connection: str
    :param force_download: force method to download
    :type force_download: bool
    :param int,list,tuple taxids: int or iterable of NCBI taxonomy identifiers (default is None = load all)
    """
    if isinstance(taxids, int):
        taxids = (taxids,)
    db = DbManager(connection)
    db.db_import_xml(urls, force_download, taxids, silent)
    db.session.close()

def set_connection(connection=defaults.sqlalchemy_connection_string_default):
    """
    Set the connection string for sqlalchemy and writes to the configuration file
    :param str connection: sqlalchemy connection string
    """
    cfp = defaults.config_file_path
    config = RawConfigParser()

    connection = connection.strip()

    if not os.path.exists(cfp):
        with open(cfp, 'w') as config_file:
            config['database'] = {'sqlalchemy_connection_string': connection}
            config.write(config_file)
            log.info('create configuration file {}'.format(cfp))
    else:
        config.read(cfp)
        config.set('database', 'sqlalchemy_connection_string', connection)
        with open(cfp, 'w') as configfile:
            config.write(configfile)

    @classmethod
    def download(cls, url=None, force_download=False):
        """Downloads uniprot_sprot.xml.gz and reldate.txt (release date information) from URL or file path

        .. note::

            only URL/path of xml.gz is needed and valid value for parameter url. URL/path for reldate.txt have to be the
            same folder
    
        :param str url: UniProt gzipped URL or file path
        :param force_download: force method to download
        :type force_download: bool
        """
        if url:
            version_url = os.path.join(os.path.dirname(url), defaults.VERSION_FILE_NAME)
        else:
            url = os.path.join(defaults.XML_DIR_NAME, defaults.SWISSPROT_FILE_NAME)
            version_url = os.path.join(defaults.XML_DIR_NAME, defaults.VERSION_FILE_NAME)
        #print(url)
        #print(version_url)
        xml_file_path = cls.get_path_to_file_from_url(url)
        version_file_path = cls.get_path_to_file_from_url(version_url)

        if force_download or not os.path.exists(xml_file_path):

            log.info('download {} and {}'.format(xml_file_path, version_file_path))

            scheme = urlsplit(url).scheme

            if scheme in ('ftp', 'http'):
                urlretrieve(version_url, version_file_path)
                #urlretrieve(url, xml_file_path)

            elif not scheme and os.path.isfile(url):
                shutil.copyfile(url, xml_file_path)
                shutil.copyfile(version_url, version_file_path)

        print(version_file_path)
        print(xml_file_path)
        return xml_file_path, version_file_path

    @classmethod
    def get_path_to_file_from_url(cls, url):
        """standard file path
        
        :param str url: download URL
        """
        file_name = urlparse(url).path.split('/')[-1]
        return os.path.join(PYUNIPROT_DATA_DIR, file_name)

set_connection('postgresql://uniport_user:uniport_user@localhost/uniprot_db')
#pyuniprot.set_mysql_connection(host='172.21.0.138', user='pyuniprot_user', passwd='PYuniprot_user#123', db='pyuniprot')
update(force_download=True)
