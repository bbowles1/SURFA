import json
import sqlite3
import pandas as pd
import logging
from importlib.resources import files
import hashlib

logger = logging.getLogger(__name__)

__all__ = ["load_schema", "validate_and_coerce", "create_table", "write_to_db"]

PANDAS_TYPE_MAP = {
    "TEXT": "object",
    "INTEGER": "int64",
    "REAL": "float64",
    "BLOB": "object",
    "NUMERIC": "float64",
}


def calculate_md5(file_path):
    """Calculate md5sum of file from path

    :param file_path: path to file
    :type file_path: str
    :return: md5 hash
    :rtype: str
    """

    # init md5 object
    md5_hash = hashlib.md5()

    with open(file_path, "rb") as f:
        # read file in 8kb chunks
        while chunk := f.read(8192):
            md5_hash.update(chunk)

    # return hexdigest
    md5sum = md5_hash.hexdigest()
    logger.debug(f"md5sum of input file {file_path} is {md5sum}.")
    return md5sum


def create_metadata_df(metadata_dict):
    """Create a pandas dataframe of build metadata

    :param metadata_dict: _description_
    :type metadata_dict: _type_
    :return: _description_
    :rtype: _type_
    """

    md_df = pd.DataFrame(list(metadata_dict.items()), columns=["input", "path"])
    md_df["md5sum"] = ""

    for name, path in metadata_dict.items():
        if path:
            md5sum = calculate_md5(path)
        else:
            logger.debug(
                f"No path was provided for input with name {name}. md5sum set to NULL."
            )
            md5sum = None

        # update df with md5
        md_df.loc[md_df.input == name, "md5sum"] = md5sum

    return md_df


def load_schema() -> dict:
    """Load schema.json from package installation dir

    :return: Schema json
    :rtype: dict
    """
    schema_path = files("surfa").joinpath("schema.json")
    logger.debug(f"Loading schema from path {schema_path}.")

    schema_text = schema_path.read_text()
    return json.loads(schema_text)


def validate_and_coerce(df: pd.DataFrame, table_schema: dict) -> pd.DataFrame:
    """Validate a DataFrame against a table schema and coerce types.
    :param df: Dataframe to convert to table.
    :type df: pd.DataFrame
    :param table_schema: dictionary of table schema.
    :type table_schema: dict
    :raises ValueError: Missing columns in input df.
    :raises ValueError: Null values inserted into non-nullable table.
    :raises TypeError: Error coercing pandas type to sql type.
    :return: Pandas dataframe with coerced types
    :rtype: pd.DataFrame
    """

    schema_columns = {col["name"] for col in table_schema["columns"]}
    df_columns = set(df.columns)

    # determine missing columns if any
    missing = schema_columns - df_columns
    if missing:
        raise ValueError(f"Table '{table_schema['name']}': missing columns {missing}")

    for col in table_schema["columns"]:
        name = col["name"]
        expected_type = PANDAS_TYPE_MAP[col["type"]]

        # check for null values if col is nullable
        if not col.get("nullable", True) and df[name].isnull().any():
            raise ValueError(
                f"Table '{table_schema['name']}': column '{name}' has nulls but is non-nullable"
            )

        # coerce type
        try:
            df[name] = df[name].astype(expected_type)
        except (ValueError, TypeError) as e:
            raise TypeError(
                f"Table '{table_schema['name']}': cannot coerce '{name}' to {expected_type}: {e}"
            )

    # Return columns in schema-defined order
    return df[[col["name"] for col in table_schema["columns"]]]


def create_table(conn: sqlite3.Connection, table_schema: dict):
    """Create a table from a schema definition

    :param conn: Existing sqlite3 connection
    :type conn: sqlite3.Connection
    :param table_schema: dictionary with format {table_name:pandas.dataFrame}
    :type table_schema: dict
    """

    columns = []
    for col in table_schema["columns"]:
        name_type_pairs = [f'"{col["name"]}"', PANDAS_TYPE_MAP[col["type"]]]
        # determine if col is primary key
        if col.get("primary_key"):
            name_type_pairs.append("PRIMARY KEY")
        # determine if col is nullable
        if not col.get("nullable", True):
            name_type_pairs.append("NOT NULL")
        columns.append(" ".join(name_type_pairs))

    if "foreign_keys" in table_schema:
        for fk in table_schema["foreign_keys"]:
            columns.append(
                f'FOREIGN KEY ("{fk["column"]}") '
                f'REFERENCES "{fk["references_table"]}" ("{fk["references_column"]}")'
            )

    col_sql = ",\n  ".join(columns)
    conn.execute(f'DROP TABLE IF EXISTS "{table_schema["name"]}"')
    conn.execute(f'CREATE TABLE "{table_schema["name"]}" (\n  {col_sql}\n)')


def write_to_db(dataframes: dict[str, pd.DataFrame], db_path: str):
    """
    Write a dict of DataFrames to SQLite using an input JSON schema.

    Args:
        dataframes: mapping of table name to DataFrame, e.g.
                    {"transcripts": transcript_df, "uorfs": uorf_df}
        schema_path: path to the JSON schema file
        db_path: path to the output SQLite database
    """

    # load schema from package src
    schema = load_schema()

    table_schemas = {t["name"]: t for t in schema["tables"]}

    # Check all provided DataFrames have a matching schema entry
    for name in dataframes.keys():
        if name not in table_schemas.keys():
            raise ValueError(f"No schema definition found for table '{name}'!")
        else:
            logger.debug(f"Found schema for {name}.")

    conn = sqlite3.connect(db_path)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        for name, df in dataframes.items():
            table_schema = table_schemas[name]
            df = validate_and_coerce(df, table_schema)
            logger.info(f"Writing table {name}.")
            create_table(conn, table_schema)
            df.to_sql(name, conn, if_exists="append", index=False)
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()
