#!/usr/bin/env zsh

cd "$(dirname "$0")"

awk '
BEGIN   {
    type["uint"]="INTEGER"
    type["int"]="INTEGER"
    type["double"]="INTEGER"
    type["str"]="TEXT"
}

/^TABLE/ 	    { print "CREATE TABLE " $2 "(" }
/^CN|^C|^LC/	{ print $2 " " type[$3] " NOT NULL," }

/^#SQL/         { print }

/^END WITH_ROWID/	    { print ");" }
/^END WITHOUT_ROWID/    { print ") WITHOUT ROWID;" }
/^END/				    { printf "\n" }

' <<< $(cat ../schema.txt | grep -v "^//") | sed 's/#SQL //g'
