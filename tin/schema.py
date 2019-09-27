# -*- coding: utf-8 -*-

"""Database schema setup."""

import logging

from psycopg2 import sql

log = logging.getLogger(__name__)


def create_schema(conn, schema):
    query_params = {'schema': schema.schema.sqlid}
    query = sql.SQL("""
    CREATE SCHEMA IF NOT EXISTS {schema}
    """).format(**query_params)
    log.debug(conn.print_query(query))
    conn.send_query(query)


def create_tables(conn, schema, epsg: int):
    query_params = {
        'tin': schema.schema + schema.table,
        'id': schema.field.id.sqlid,
        'geom': schema.field.geom.sqlid,
        'star': schema.field.star.sqlid,
        'epsg': sql.Literal(epsg)
    }
    query = sql.SQL("""
    CREATE TABLE IF NOT EXISTS {tin} ( 
        {id} bigint, 
        {geom} geometry(POINTZ, {epsg}), 
        {star} bigint[]
    )
    """).format(**query_params)
    log.debug(conn.print_query(query))
    conn.send_query(query)


def create_relations(conn, schema, epsg):
    conn.check_postgis()
    create_schema(conn, schema)
    create_tables(conn, schema, epsg)


def drop_relations(conn, schema):
    query = sql.SQL("""
    DROP SCHEMA IF EXISTS {} CASCADE
    """).format(schema.schema.sqlid)
    log.debug(conn.print_query(query))
    conn.send_query(query)