all: ./sql/sql-structs.hh ./sql/sql-str-query.hpp ./cbs.hpp

./sql/sql-structs.hh: ./schema.txt ./sql/sql-structs.sh
	./sql/sql-structs.sh 

./sql/sql-str-query.hpp: ./schema.txt ./sql/sql-str-query.sh 
	./sql/sql-str-query.sh 

./cbs.hpp : ./schema.txt ./cbs.hpp.sh
	./cbs.hpp.sh
