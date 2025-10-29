import zipfile
import sqlite3

from mddb_workflow.utils.auxiliar import InputError, save_json
from mddb_workflow.utils.constants import OUTPUT_PROVENANCE_FILENAME
from mddb_workflow.utils.type_hints import *

def produce_provenance (
    output_directory : str,
    aiida_data_file : Optional['File'],
):
    """Produce a provenance file containing AiiDA data adapted for our database"""

    # The AiiDA exported file is a zip-compressed sqlite file
    if not aiida_data_file or not aiida_data_file.exists:
        print(' There is no AiiDA data')
        return

    # So the first step is extracting all zip contents
    try:
        with zipfile.ZipFile(aiida_data_file.path, 'r') as zip_ref:
            zip_ref.extractall(output_directory)
    except zipfile.BadZipFile:
        raise InputError('AiiDA data file must be a zip-compressed sqlite.\n' +
            '  The file you provided is not zip-compressed: ' + aiida_data_file.path)
    except:
        raise RuntimeError(f'Something went wrong while decompressing zip file {aiida_data_file.path}')

    # Now find the db file among decompressed files and parse it
    # DANI: Asumo que siempre se llamará igual, pero quien sabe
    db_filepath = f'{output_directory}/db.sqlite3'

    # Setup the database
    connection = sqlite3.connect(db_filepath)
    cursor = connection.cursor()

    # Mine the target tables
    # DANI: El día que hice esto el resto de tablas estaban vacías
    target_tables = { 'db_dbcomputer', 'db_dbuser', 'db_dbnode', 'db_dblink' }
    tables = {}
    for table_name in target_tables:
        cursor.execute(f'select * from {table_name}')
        headers = [desc[0] for desc in cursor.description]
        all_rows = cursor.fetchall()
        tables[table_name] = { 'headers': headers, 'rows': all_rows }

    # Save data in JSON format
    output_filepath = f'{output_directory}/{OUTPUT_PROVENANCE_FILENAME}'
    save_json(tables, output_filepath)