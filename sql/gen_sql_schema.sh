#!/usr/bin/env zsh

typeset -A types
types[uint]="INTEGER"
types[int]="INTEGER"
types[double]="INTEGER"
types[str]="TEXT"

while read line; do
    case $line in
        TABLE*)
            table_name=$(echo $line | cut -d" " -f2)
            echo "create table $table_name("
            ;;
        "#END WITH_ROWID")
            echo ");"
            ;;
        "#END WITHOUT_ROWID")
            echo ") WITHOUT ROWID;"
            ;;
        \#SQL*)
            echo $line | cut -d" " -f2-
            ;;
        *" "*)
            var_name=$(echo $line | cut -d" " -f1)
            var_type=$(echo $line | cut -d" " -f2)
            echo $var_name ${types[$var_type]} NOT NULL,
            ;;
        "")
            echo
            ;;
    esac
done < $1
