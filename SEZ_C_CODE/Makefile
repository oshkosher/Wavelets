default: main_comp main_ucomp

CC = gcc -g

main_comp: sez_alloc.h sez_wcomp.h zlib_wrapper.h \
  main_comp.c sez_alloc.c sym_fwt_2d.c sez_wcomp.c zlib_wrapper.c comp.c
	$(CC) main_comp.c sez_alloc.c sym_fwt_2d.c sez_wcomp.c \
	  zlib_wrapper.c comp.c -lz -lm -o $@


main_ucomp: sez_alloc.h sez_wcomp.h zlib_wrapper.h \
  main_ucomp.c sez_alloc.c sym_iwt_2d.c sez_wcomp.c zlib_wrapper.c comp.c
	$(CC) main_ucomp.c sez_alloc.c sym_iwt_2d.c sez_wcomp.c \
	  zlib_wrapper.c ucomp.c -lz -lm -o $@

clean:
	rm -f main_comp main_ucomp
