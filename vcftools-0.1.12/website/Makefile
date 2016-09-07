
FILES=documentation downloads examples htslib index license links man_latest perl_examples perl_module specs

install:
	    @for i in $(FILES); do echo "php src/index.php $$i > $$i.html"; php src/index.php $$i > $$i.html; done;
