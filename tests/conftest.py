# -*- coding: utf-8 -*-

"""pytest configuration"""


import os
from pathlib import Path
import pytest
from tin import db, main

#-------------------------------------------------------------------- testing DB
@pytest.fixture('session')
def cfg():
    config = '/home/balazs/Development/TIN_scale_up/tin_config.yml'
    with open(config, 'r') as fo:
        c = main.parse_config(fo)
        yield c


@pytest.fixture('session')
def tin_db(cfg):
    conn = db.Db(**cfg['database'])
    yield conn
    conn.close()



@pytest.fixture('session')
def tin_schema(cfg):
    yield db.Schema(cfg['tin_schema'])


#-------------------------------------------------------------------- directory
@pytest.fixture('session')
def t_dir():
    """tests directory"""
    yield Path(__file__).parent


@pytest.fixture('session')
def data_dir(t_dir):
    yield t_dir / 'data'


@pytest.fixture('session')
def root_dir(t_dir):
    yield t_dir.parent


@pytest.fixture('session')
def package_dir(root_dir):
    yield root_dir / 'tin'
