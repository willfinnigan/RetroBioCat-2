

class DB_Query_SQLite():

    def __init__(self, database):
        self.database = database
        self.conn = self.database.conn

    def smi_exists(self, smiles, table):
        cursor = self.conn.cursor()
        query = f"SELECT EXISTS(SELECT 1 FROM {table} WHERE SMILES='{smiles}');"
        cursor.execute(query)
        result = cursor.fetchall()[0][0]
        cursor.close()
        return result

    def vendor_rows(self, vendor, table):
        cursor = self.conn.cursor()
        cmd = f"select count(*) from {table} where {vendor}_id is not NULL"
        cursor.execute(cmd)
        result = cursor.fetchone()
        cursor.close()

        return result

    def get_vendor_counts(self, table, vendors):
        for v in vendors:
            num = self.vendor_rows(v, table)

    def smiles_lookup(self, smiles, table, vendors=None):
        if vendors == []:
            return None

        cursor = self.conn.cursor()
        cmd = f"select * from {table} where smiles='{smiles}'"

        if vendors is not None:
            cmd += " and ("
            for i, vendor in enumerate(vendors):
                cmd += f"{vendor}_id is not NULL"

                if i != len(vendors) - 1:
                    cmd += " or "
                else:
                    cmd += ")"

        cursor.execute(cmd)
        result = cursor.fetchone()
        cursor.close()

        return result

    def get_column_names(self, table):
        cursor = self.conn.cursor()
        cmd = f"select * from pragma_table_info('{table}')"
        cursor.execute(cmd)
        result = cursor.fetchall()
        names = [i[1] for i in result]
        cursor.close()
        return names

