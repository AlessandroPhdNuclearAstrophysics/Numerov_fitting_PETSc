.PHONY = all clean

RM := rm -vf
MAKE := make
TARGET := fit_potential.x

OUTPUT_DIR := output

all: 
	@$(MAKE) -j -C build all

run: all
	./bin/$(TARGET) -tao_monitor_short -tao_max_it 10000 -tao_type pounders -tao_gatol 1.e-8

compare:
	@echo  "@with g0 \n\
	@		  title \"Compare fitted solution to AV18 solution for 1P1\" \n\
	@		  xaxis  label \"k\\S2\\N (fm\\S-2\\N)\" \n\
	@		  yaxis  label \"k\\S3\\N cot delta (fm\\S-3\\N)\" \n\
	@     s0 legend \"AV18\" \n\
	@     s1 legend \"fitted EFT-pless\" \n\
	@     legend on \n\
	@		  legend 0.8, 0.8 \n\
	@     s0 line linewidth 2 \n\
	@     s1 line linewidth 2 \n\
	@     s1 line linestyle 3 \n\
	"> config.agr
	@sed -i 's/delta/\\xd\\f{}/g' config.agr
	xmgrace $(OUTPUT_DIR)/*.dat config.agr -saveall $(OUTPUT_DIR)/compare.agr
	@$(RM) -v config.agr


clean:
	@$(MAKE) -C build clean

del_out:
	@$(RM) -v output/*.dat
