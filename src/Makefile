LIBNAME=rmutil
LD=ld

OBJS=	dist.o gettvc.o inthaz.o romberg.o toms614.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(LD) $(SHLIBLDFLAGS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
