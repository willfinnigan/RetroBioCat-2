import pandas as pd
import csv

from rbc2.reaction_evaluation.sqlite_source_mol.connect_sqlitedb import SQLite_Database


class DB_Creator_SQLite():

    def __init__(self, database):
        self.database = database
        self.conn = self.database.conn
        self.engine = self.database.sqlalchemy_engine()

    def load_from_csv(self, csv_path, tablename, columns=None, dtypes=None):

        if columns is None:
            with open(csv_path, 'r') as infile:
                reader = csv.DictReader(infile)
                columns = reader.fieldnames[1:]

        if dtypes is None:  # if no dtypes set, then everything is a string
            dtypes = {k: str for k in columns}

        df = pd.read_csv(csv_path, index_col=0, dtype=dtypes)

        if 'SMILES' in df.columns:  # make SMILES lowercase if necessary
            df.rename(columns={"SMILES": "smiles"}, inplace=True)

        self._data_to_table(df, tablename)
        self._create_index(tablename, 'smiles')

    def _data_to_table(self, df, table):
        print("data to sqlite..")
        df.to_sql(table, con=self.engine, if_exists='replace')
        print('done')

    def _create_index(self, table, col_for_index):
        cursor = self.conn.cursor()
        cmd = f"CREATE INDEX {table}_idx_{col_for_index} ON {table} ({col_for_index});"
        cursor.execute(cmd)
        self.conn.commit()
        cursor.close()
        print('done')


if __name__ == '__main__':
    db = SQLite_Database('source_mols.db')
    db_creator = DB_Creator_SQLite(db)
    db_creator.load_from_csv('new_building_blocks_final.csv', 'building_blocks')
    db_creator.load_from_csv('metabolites_final.csv', 'metabolites')

