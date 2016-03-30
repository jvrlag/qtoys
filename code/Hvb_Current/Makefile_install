# HVB makefile

## Include configurations
confs=make.inc
ifeq ($(wildcard $(confs)),$(confs))
include $(confs)
else
$(error File $(confs) required)
endif
##


all:
	@echo "Use:"
	@echo "make test"
	@echo "make install [may require root privileges]"


.SILENT:test
test:
	@echo "TESTING"
	@echo -ne "\nSource dir: $(SRCD) : "
	test -d $(SRCD) && echo -n "found : " || echo -n "NOT FOUND, "
	test -w $(SRCD) && echo "writeable"|| echo "NON-WRITEABLE"
	@echo -n "Directory lib: $(LIBD) : "
	test -d $(LIBD) && echo -n "found : " || echo -n "NOT FOUND, "
	test -w $(LIBD) && echo "writeable"|| echo "NON-WRITEABLE"
	@echo -n "Directory include: $(INCD) : "
	test -d $(INCD) && echo -n "found : " || echo -n "NOT FOUND, "
	test -w $(INCD) && echo "writeable"|| echo "NON-WRITEABLE"
	@echo -ne "\nDirectory : $(INCD)$(NAME)$(VERSION) : "
	test -d $(INCD)$(NAME)$(VERSION) && echo "Exists: CHECK" || echo "to be created"
	@echo -ne "External libreries\n"
	@echo -n "Imlib2 : "
	test -e /usr/include/Imlib2.h && echo "OK" || echo "NOT FOUND"
	@echo -n "X11/Xlib : "
	test -e /usr/include/X11/Xlib.h && echo "OK" || echo "NOT FOUND"
	@echo -n "X11/Xutil : "
	test -e /usr/include/X11/Xutil.h && echo "OK" || echo "NOT FOUND"
	@echo -n "X11/keysym.h : "
	test -e /usr/include/X11/keysym.h && echo "OK" || echo "NOT FOUND"
	@echo -n "Soft-Link : $(INCD)$(NAME) : "
	test -e $(INCD)$(NAME) && echo "Exists: CHECK" || echo "to be created"
	@echo -n "Library : $(LIBD)$(LNAME) : "
	test -e $(LIBD)$(LNAME) && echo "Exists: CHECK" || echo "to be created"
	@echo -n "Soft-Link : $(LIBD)$(GNAME) : "
	test -e $(LIBD)$(GNAME) && echo "Exists: CHECK" || echo "to be created"
	@echo -ne "\nLibrary-files to be included: "
	@echo $(HVBH)
	@echo -ne "\nSource-files to be compiled: "
	@echo $(HVBSRC)
	@echo -ne "\nObjects to be generated: "
	@echo $(HVBS)


install: objects
	@echo "INSTALLING"
	@mkdir -p $(INCD)$(NAME)$(VERSION)
	@$(foreach i,$(HVBH), cp -vi $(SRCD)$(i) $(INCD)$(NAME)$(VERSION);)
	@echo "BUILDING LIBRARY"
	ar rcs $(LIBD)$(LNAME) $(HVBS)
	@test -e $(INCD)$(NAME) && rm $(INCD)$(NAME) || echo "$(NAME)"
	cd $(INCD) && ln -s $(NAME)$(VERSION) $(NAME)
	@test -e $(LIBD)$(GNAME) && rm $(LIBD)$(GNAME) || echo "$(GNAME)"
	cd $(LIBD) && ln -s $(LNAME) $(GNAME)


objects: $(HVBS)
%.o: %.cc
	$(CC) -c $< $(FLAGS) -o $@


.PHONY:clean
clean:
	$(foreach i,$(wildcard *.o), rm $(i);)
	$(foreach i,$(wildcard *.a), rm $(i);)
