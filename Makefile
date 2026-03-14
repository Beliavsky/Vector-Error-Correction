executables = xvecm_gfort.exe xvecm_rank_gfort.exe xvecm_rank_gap_gfort.exe
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private
obj    = kind.o random.o linear.o vecm.o xvecm.o xvecm_rank.o xvecm_rank_gap.o

all: $(executables)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

xvecm_gfort.exe: kind.o random.o linear.o vecm.o xvecm.o
	$(FC) -o xvecm_gfort.exe kind.o random.o linear.o vecm.o xvecm.o $(FFLAGS)

xvecm_rank_gfort.exe: kind.o random.o linear.o vecm.o xvecm_rank.o
	$(FC) -o xvecm_rank_gfort.exe kind.o random.o linear.o vecm.o xvecm_rank.o $(FFLAGS)

xvecm_rank_gap_gfort.exe: kind.o random.o linear.o vecm.o xvecm_rank_gap.o
	$(FC) -o xvecm_rank_gap_gfort.exe kind.o random.o linear.o vecm.o xvecm_rank_gap.o $(FFLAGS)

run: $(executables)
	./xvecm_gfort.exe
	./xvecm_rank_gfort.exe
	./xvecm_rank_gap_gfort.exe

clean:
	rm -f $(executables) $(obj)

