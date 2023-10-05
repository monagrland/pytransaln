help:
	@echo "make build        Install package with setuptools and run unit tests"
	@echo "make benchmark    Download test data (from transAlign benchmark set) and perform alignments"
	@echo "make clean        Delete benchmark test run output"

build: pyproject.toml
	pip install .
	python -m unittest -v pytransaln.testmodule

clean:
	rm benchmark/test.*

benchmark: benchmark/test.cons.RAG2.aln benchmark/test.cons.RBP3.aln benchmark/test.each.RAG2.aln benchmark/test.each.RBP3.aln benchmark/test.each.BDNF.aln benchmark/test.each_notermstop.BDNF.aln

benchmark/test.cons.%.aln: benchmark/benchMark_data/%_unaligned.fas
	pytransaln --input $< --code 1 align --how cons --threads 8 --out_aa $@.aa --out_bad $@.bad --out_aln_aa $@.aln_aa --out_aln_nt $@.aln_nt --out_aln_nt_aug $@ --out_bad_fs_report $@.fs_report

benchmark/test.each_notermstop.%.aln: benchmark/benchMark_data/%_unaligned.fas
	pytransaln --input $< --code 1 --ignore_terminal_stop align --how each --threads 8 --out_aa $@.aa --out_bad $@.bad --out_aln_aa $@.aln_aa --out_aln_nt $@.aln_nt --out_aln_nt_aug $@ --out_bad_fs_report $@.fs_report

benchmark/test.each.%.aln: benchmark/benchMark_data/%_unaligned.fas
	pytransaln --input $< --code 1 align --how each --threads 8 --out_aa $@.aa --out_bad $@.bad --out_aln_aa $@.aln_aa --out_aln_nt $@.aln_nt --out_aln_nt_aug $@ --out_bad_fs_report $@.fs_report

benchmark/benchMark_data/%_unaligned.fas: benchmark/benchMark_data.tar.gz
	tar -C ./$(<D) -xzf $< benchMark_data/$(@F)

benchmark/benchMark_data.tar.gz:
	mkdir -p benchmark;
	wget --no-clobber -P ./$(@D) https://uol.de/f/5/inst/biologie/ag/systematik/download/Programs/benchMark_data.tar.gz    

