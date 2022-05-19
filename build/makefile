# author: mhnrchs

# compiler
CC = g++
# compiler flags
CFLAGS = -std=c++14 -Wall -g -O3

all: full filter_lookup

# executables
full: full.o vcf2binmat.o sel_coeff_est.o create_hist.o utils.o pbwt.o filter_blocks.o
	$(CC) $(CFLAGS) full.o vcf2binmat.o sel_coeff_est.o create_hist.o utils.o pbwt.o filter_blocks.o -o full

# out files
full.o: full.cpp
	$(CC) $(CFLAGS) -c full.cpp
vcf2binmat.o: vcf2binmat.cpp vcf2binmat.h
	$(CC) $(CFLAGS) -c vcf2binmat.cpp
sel_coeff_est.o: sel_coeff_est.cpp
	$(CC) $(CFLAGS) -c sel_coeff_est.cpp
create_hist.o: create_hist.cpp create_hist.h
	$(CC) $(CFLAGS) -c create_hist.cpp
pbwt.o: pbwt.cpp pbwt.h
	$(CC) $(CFLAGS) -c pbwt.cpp
utils.o: utils.cpp utils.h
	$(CC) $(CFLAGS) -c utils.cpp
filter_blocks.o: filter_blocks.cpp filter_blocks.h
	$(CC) $(CFLAGS) -c filter_blocks.cpp

filter_lookup: stat_filter/filter_lookup.o stat_filter/brent.o utils.o
	$(CC) $(CFLAGS) stat_filter/filter_lookup.o utils.o stat_filter/brent.o -o filter_lookup
stat_filter/filter_lookup.o: stat_filter/filter_lookup.cpp
	$(CC) $(CFLAGS) -c stat_filter/filter_lookup.cpp -o stat_filter/filter_lookup.o
stat_filter/brent.o: stat_filter/brent.cpp stat_filter/brent.hpp
	$(CC) $(CFLAGS) -c stat_filter/brent.cpp -o stat_filter/brent.o

# standalones
standalones: vcf2binmat pbwt sel_coeff_est filter_blocks create_hist

vcf2binmat: standalones/vcf2binmat_standalone.o vcf2binmat.o utils.o
	$(CC) $(CFLAGS) standalones/vcf2binmat_standalone.o vcf2binmat.o utils.o -o vcf2binmat
standalones/vcf2binmat_standalone.o: standalones/vcf2binmat_standalone.cpp
	$(CC) $(CFLAGS) -c standalones/vcf2binmat_standalone.cpp -o standalones/vcf2binmat_standalone.o

pbwt: standalones/pbwt_standalone.o pbwt.o utils.o
	$(CC) $(CFLAGS) standalones/pbwt_standalone.o pbwt.o utils.o -o pbwt
standalones/pbwt_standalone.o: standalones/pbwt_standalone.cpp
	$(CC) $(CFLAGS) -c standalones/pbwt_standalone.cpp -o standalones/pbwt_standalone.o

sel_coeff_est: standalones/sel_coeff_est_standalone.o sel_coeff_est.o utils.o
	$(CC) $(CFLAGS) standalones/sel_coeff_est_standalone.o sel_coeff_est.o utils.o -o sel_coeff_est
standalones/sel_coeff_est_standalone.o: standalones/sel_coeff_est_standalone.cpp
	$(CC) $(CFLAGS) -c standalones/sel_coeff_est_standalone.cpp -o standalones/sel_coeff_est_standalone.o

filter_blocks: standalones/filter_blocks_standalone.o filter_blocks.o utils.o
	$(CC) $(CFLAGS) standalones/filter_blocks_standalone.o filter_blocks.o utils.o -o filter_blocks
standalones/filter_blocks_standalone.o: standalones/filter_blocks_standalone.cpp
	$(CC) $(CFLAGS) -c standalones/filter_blocks_standalone.cpp -o standalones/filter_blocks_standalone.o

create_hist: standalones/create_hist_standalone.o create_hist.o utils.o
	$(CC) $(CFLAGS) standalones/create_hist_standalone.o create_hist.o utils.o -o create_hist
standalones/create_hist_standalone.o: standalones/create_hist_standalone.cpp
	$(CC) $(CFLAGS) -c standalones/create_hist_standalone.cpp -o standalones/create_hist_standalone.o

clean:
	rm *.o
	rm stat_filter/*.o
	rm standalones/*.o
