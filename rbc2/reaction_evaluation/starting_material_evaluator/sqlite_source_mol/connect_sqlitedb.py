import sqlite3
import sqlalchemy
from rbc2.configs.data_path import path_to_data_folder

data_folder = f'{path_to_data_folder}/buyability'
DEFAULT_DB_PATH = f'{data_folder}/source_mols.db'

class SQLite_Database():

    def __init__(self, path):
        if path is None:
            path = DEFAULT_DB_PATH
            print(path)

        self.engine = None
        self.path = path
        self.conn = self.connect()

    def connect(self):
        return sqlite3.connect(self.path)

    def sqlalchemy_engine(self):
        if self.engine is None:
            conn_string = f"sqlite:///{self.path}"
            self.engine = sqlalchemy.create_engine(conn_string)
        return self.engine