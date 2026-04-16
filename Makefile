clean:
	find . -type f \( -name "*.log" -o -name "*.aux" -o -name "*.synctex.gz" -o -name "*.fdb_latexmk" -o -name "*.fls" -o -name "*.out" \) -delete

commit: 
	bash autocommit.sh