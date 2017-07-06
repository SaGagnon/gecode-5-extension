#!/usr/bin/env zsh

cd $(dirname "$0")

function replace_between {
    awk -v TO_INSERT="$1" -v BEGIN_="$2" -v END_="$3" '
    BEGIN { p=1 }
    $0~BEGIN_ { print; print TO_INSERT; p=0}
    $0~END_ { p=1 }
    p' < cbs.hpp | sponge cbs.hpp
}

TXT=$(cat schema.txt | sed -n '/C max_sd/,/LC /p')

################################################################################


TRANS=$(awk '
BEGIN {
 type["double"]="double"
 type["uint"]="unsigned int"
 type["int"]="int"
}
{ print "static " type[$3] " " $2 "[SIZE];" }
' <<< $TXT)

BEGIN="//#BEGIN_ATTRIBUTES"
END="//#END_ATTRIBUTES"

replace_between $TRANS $BEGIN $END
