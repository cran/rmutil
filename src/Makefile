LIBNAME=rmutil

OBJS=	dist.o gettvc.o inthaz.o romberg.o toms614.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(SHLIBLD) $(SHLIBLDFLAGS) $(FLIBS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
