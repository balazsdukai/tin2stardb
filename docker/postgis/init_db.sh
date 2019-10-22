#!/bin/bash
set -e

psql -p 5432 -h localhost -U postgres -d postgres -c "create role tin_tester with superuser login password 'tin_tester';"
psql -p 5432 -h localhost -U postgres -d postgres -c "create database tin_db with owner tin_tester;"
