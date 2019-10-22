#!/bin/bash
set -e

psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" <<-EOSQL
    CREATE USER tin_tester WITH LOGIN PASSWORD 'tin_tester_pw';
    CREATE DATABASE tin_db;
    GRANT ALL PRIVILEGES ON DATABASE tin_db TO tin_tester;
EOSQL