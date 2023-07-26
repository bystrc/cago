# Start of the makefile
# Defining variables
objects = pdbx_dsp.o cagov2.o
f95comp = mpif90
# Makefile
xcago: $(objects)
	$(f95comp) -o xcago $(objects)
pdbx_dsp.mod: pdbx_dsp.o pdbx_dsp.f95
	$(f95comp) -c pdbx_dsp.f95
pdbx_dsp.o: pdbx_dsp.f95
	$(f95comp) -c pdbx_dsp.f95
cagov2.o: pdbx_dsp.mod cagov2.f95
	$(f95comp) -c cagov2.f95

# Cleaning everything
clean:
	rm pdbx_dsp.mod xcago
	rm $(objects)
# End of the makefile
